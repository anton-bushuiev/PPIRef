import os
import subprocess
import shutil
import shlex
from typing import Union, Literal, Optional
from pathlib import Path

import pandas as pd

from ppiref.utils.residue import Residue
from ppiref.utils.pdb import get_n_models, get_first_model
from ppiref.utils.misc import random_id
from ppiref.definitions import DR_SASA_PATH


class DR_SASA:
    def __init__(
        self,
        path: Union[Path, str] = DR_SASA_PATH,
        tmp_dir: Union[Path, str] = None,
        verbose: bool = False,
        auto_clean: Literal['lazy', 'instant', 'none'] = 'lazy'
    ) -> None:
        """Python wrapper for the dr_sasa software (https://github.com/nioroso-x3/dr_sasa_n).
        
        Args:
            path: Path to dr_sasa executable. Defaults to DR_SASA_PATH.
            tmp_dir: Path to temporary dir to store outputs. Defaults to None.
        """
        self.path = Path(path)
        self.tmp_dir = Path(tmp_dir) if tmp_dir else Path(f'./.dr_sasa_tmp_dir_{random_id(20)}')
        self.tmp_dir = self.tmp_dir.resolve()  # Make absolute
        self.tmp_dir.mkdir(exist_ok=False, parents=True)
        self.verbose = verbose
        self.auto_clean = auto_clean

    def __del__(self):
        """Clean on destruction.
        """
        if self.auto_clean == 'lazy':
            self.clean()

    def __call__(
        self,
        pdb_path: Union[Path, str],
        partners: tuple[str]
    ) -> tuple[float, set]:
        """Call dr_sasa executable and return parsed output.

        TODO Extend to classify buried residues according to Levy 2010

        Args:
            pdb_path: Path to input .pdb file
            partners: _description_

        Returns:
            _description_
        """
        # Prepare args
        pdb_path = Path(pdb_path)
        pdb_stem = pdb_path.stem
        a, b = partners
        remove_pdb = False

        # Extract first model only
        if get_n_models(pdb_path) > 1:
            pdb_path_first_model = self.tmp_dir / f'{pdb_stem}.pdb'
            get_first_model(pdb_path, pdb_path_first_model)
            pdb_path = pdb_path_first_model
            remove_pdb = True

        # Execute dr_sasa from temporary directory
        command = f'{self.path} -m 1 -i {pdb_path} -chain {a} -chain {b}'
        command = shlex.split(command)
        subprocess.run(command, cwd=self.tmp_dir, check=True, capture_output=not self.verbose)

        # Read pairwise delta SASA matrix
        matrix_path = self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.matrix.{a}{b}.by_res.tsv'
        df_pairwise_dsasa = pd.read_csv(matrix_path, sep='\t', index_col=0)

        # Parsed output
        # Get Buried surface area (BSA)
        bsa = df_pairwise_dsasa.sum().sum() / 2
        # Get buried residues (dSASA > 0)
        df_nonzero_dsasa = df_pairwise_dsasa.stack().reset_index()
        buried_residues = set(df_nonzero_dsasa.iloc[:, :2].to_numpy().flatten())
        buried_residues = set(map(self.parse_residue, buried_residues))

        # Clean
        if self.auto_clean == 'instant':
            self.clean(pdb_stem, partners)
        if remove_pdb:
            pdb_path.unlink()

        return buried_residues, bsa

    def clean(self, pdb_stem: Optional[str] = None, partners: Optional[tuple[str]] = None) -> None:
        """Clean whole temporary directory or one-run files.
        """
        if hasattr(self, 'tmp_dir') and self.tmp_dir.is_dir():
            if pdb_stem is None and partners is None:  # Delete all
                shutil.rmtree(self.tmp_dir)
            else:  # Delete one-run files
                a, b = partners
                files = [
                    self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.by_atom.tsv',
                    self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.by_res.tsv',
                    self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.datmasa',
                    self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.matrix.{a}{b}.by_atom.tsv',
                    self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.matrix.{a}{b}.by_res.tsv',
                    self.tmp_dir / f'{pdb_stem}.{a}_vs_{b}.overlaps',
                    self.tmp_dir / f'{pdb_stem}.{b}_vs_{a}.by_atom.tsv',
                    self.tmp_dir / f'{pdb_stem}.{b}_vs_{a}.by_res.tsv',
                    self.tmp_dir / f'{pdb_stem}.asa.pdb',
                    self.tmp_dir / f'{pdb_stem}.atmasa',
                    self.tmp_dir / f'{pdb_stem}.dsasa.pdb'
                ]
                for file in files:
                    file.unlink(missing_ok=True)

    @staticmethod
    def parse_residue(res: str) -> Residue:
        """Parse residue in a dr_sasa format into a Residue namedtuple.

        E.g. 'VAL/M/14A' -> Residue(chain_id='M', residue_number=14, insertion='A').
        """
        _, chain, pos = res.split('/')
        num = ''.join([i for i in pos if i.isdigit() or i == '-'])
        ins = pos[len(num):]
        if not pos.startswith(num):
            raise ValueError(f'Corrupted input residue {res}.')
        return Residue(chain, int(num), ins)
