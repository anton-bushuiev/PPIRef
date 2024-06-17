import os
import copy
import itertools
import traceback
import datetime
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Union, Literal, Sequence, Iterable, Optional
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn.metrics
import scipy.linalg
from tqdm import tqdm
from Bio.PDB import parse_pdb_header
from biopandas.pdb import PandasPdb

from ppiref.surface import DR_SASA
from ppiref.utils.residue import Residue, BASE_AMINO_ACIDS_3
from ppiref.utils.pdb import get_first_model
from ppiref.utils.ppipath import NOPPI_EXTENSION, path_to_pdb_id
from ppiref.utils.misc import list_to_chunks, get_partition, random_id
from ppiref.definitions import PPIREF_NAME, PPIREF_URL


INTERFACE_KIND_TYPE = Literal['heavy', 'bsa']


class PPIExtractor:
    def __init__(
        self,
        out_dir: Union[Path, str],
        kind: INTERFACE_KIND_TYPE = 'heavy',
        radius: float = 10.,
        expansion_radius: float = 0.,
        bsa: bool = False,
        join: bool = False,
        nest_out_dir: bool = True,
        max_workers: int = os.cpu_count() - 2,
        chunk_size: int = 1,
        verbose: bool = False,
        noppi_files: bool = True,
        input_format: Literal['pdb', 'haddock'] = 'pdb'
    ) -> None:
        """_summary_

        Note: With kind=='bsa' and pre-selection with radius ... TODO

        Args:
            out_dir: _description_
            kind: _description_. Defaults to 'heavy'.
            radius: _description_. Defaults to 10..
            bsa: Calculate and write BSA to output .pdb files. Defaults to False because 
                computationally expensive.
            join: TODO ... . Enables to extract beyond dimeric interfaces.
        """
        self.out_dir = Path(out_dir)
        self.kind = kind 
        self.radius = radius
        self.expansion_radius = expansion_radius
        self.bsa = bsa
        self.join = join
        self.nest_out_dir = nest_out_dir
        self.max_workers = max_workers
        self.chunk_size = chunk_size
        self.verbose = verbose
        self.noppi_files = noppi_files
        self.input_format = input_format

        if self._requires_bsa():
            self.dr_sasa = DR_SASA()

    def extract(
        self,
        pdb_path: Union[Path, str],
        partners: Optional[Iterable[str]] = None
    ) -> None:
        """Extract interfaces from `pdb_path` and write as .pdb files.

        Args:
            pdb_path: Path to .pdb file to extract PPIs from.
            partners: If not None the interfaces are extracted only between specified chains.
                In combination with self.join all interfaces between chains in `partners` are
                extracted and written into a single file. If partners is None all dimeric interfaces
                are extracted. Defaults to None.
        """
        # Preprocess args
        pdb_path = Path(pdb_path)
        pdb_id = path_to_pdb_id(pdb_path)
        pdb_stem = pdb_path.stem
        out_dir = self.out_dir / pdb_id[1:3] if self.nest_out_dir else self.out_dir
        out_dir.mkdir(exist_ok=True, parents=True)

        # Parse PDB header
        header = parse_pdb_header(pdb_path)
        pdb_stats = {
            'RESOLUTION': str(header['resolution']) + ' A',
            'STRUCTURE METHOD': header['structure_method'],
            'DEPOSITION DATE': header.get('deposition_date', None),
            'RELEASE DATE': header.get('release_date', None)
        }

        # Get first model and preprocess atom data frame
        ppdb_df = get_first_model(pdb_path)
        ppdb_df.df['ATOM'] = self._preprocess_atom_df(ppdb_df.df['ATOM'], chains=partners)

        # Get atom data frame depending on extractor type
        if self.kind in ('heavy', 'bsa'):
            atom_df = ppdb_df.get(s='heavy', records=('ATOM',))
        else:
            raise ValueError(f'Unknown `self.kind` value \'{self.kind}\'.')
        res_ids = set(atom_df['res_id'])

        # Return if less than two chains in file
        chains = atom_df['chain_id'].unique()
        if partners is not None:
            assert set(partners).issubset(chains), f'Partners {set(partners) - set(chains)} not in PDB file.'
        if len(chains) < 2:
            self._no_ppi_exit(out_dir, pdb_stem)
            return

        # Obtain inter-chain residue-residue adjacency matrix
        radius_adj = self._get_radius_adjacency(atom_df, self.radius)
        intra_adj = self._get_intra_adjacency(atom_df)
        radius_adj &= ~intra_adj

        # Extract all dimeric interfaces
        interfaces = {}
        for chains in itertools.combinations(chains, 2):
            a, b = chains
            interface = self._get_interacting_residues(radius_adj, a, b)
            if len(interface):
                interfaces[(a, b)] = interface

        # Calculate buried surface of extracted raidus-based interfaces
        bsas = {}
        if self.kind == 'bsa' or self.bsa:
            for ipartners in list(interfaces.keys()):
                buried_residues, bsa = self.dr_sasa(pdb_path, ipartners)
                bsas[ipartners] = bsa
                if self.kind == 'bsa':
                    if bsa > 0:
                        buried_residues &= res_ids  # intersect with valid residues
                        interfaces[ipartners] = buried_residues
                    else:
                        del interfaces[ipartners]

        # Expand interfaces
        expansion_adj = self._get_radius_adjacency(atom_df, self.expansion_radius)
        expansion_adj &= intra_adj
        for ipartners, interface in interfaces.items():
            interface_expansion_adj = expansion_adj.loc[pd.Series(list(interface)).to_numpy(), :]
            for p in ipartners:
                interface.update(self._get_interacting_residues(interface_expansion_adj, p, p))

        # Join interfaces
        if self.join:
            if partners is not None:
                joined_interface = set().union(*list(interfaces.values()))
                interfaces = {tuple(partners): joined_interface}
                bsas = {tuple(partners): sum(bsas.values())}
            else:
                # TODO Automatic joining based on shared residues when partners are not specified
                raise NotImplementedError()

        # Write .<pdb_id>.noppi file if no interfaces
        if not len(interfaces):
            self._no_ppi_exit(out_dir, pdb_stem)
            return

        # Write interfaces to .pdb files
        for ipartners in interfaces.keys():
            interface = interfaces[ipartners]
            interface_stats = copy.copy(pdb_stats)
            if self.bsa:
                interface_stats['BSA'] = f'{bsas[ipartners]} A^2'
            suffix = '_'.join(ipartners)
            self._write_to_pdb(
                interface,
                atom_df,
                out_dir / f'{pdb_stem}_{suffix}.pdb',
                stats=interface_stats
            )

    def _no_ppi_exit(self, out_dir: Path, pdb_stem: str):
        if self.noppi_files:
            noppi_file = out_dir / f'.{pdb_stem}{NOPPI_EXTENSION}'
            noppi_file.touch()

    def extract_parallel(
        self,
        in_dir: Union[Path, str],
        in_file_pattern: Optional[str] = '.*\.pdb$',
        partition: tuple[float] = (0., 1.),
        resume: bool = True
    ) -> None:
        """Extract interafces from all .pdb files in `in_dir` in parallel.

        Args:
            in_dir: Input directory with .pdb files to extract PPIs from.
            in_file_pattern: Pattern to match all input files in `in_dir`. Defaults to '*.pdb'.
            partition: Fractional partition of input files to process. For example, (0., 0.5) will
                process the first half of the files. This isuseful, when extracting in the data 
                parallel way across multiple nodes. Defaults to (0., 1.).
            resume: If set to True, will check what .pdb files were already processed based on the
                resulting files in the output directory, and skip them for processing. Defaults to
                True.
        """
        in_dir = Path(in_dir)

        # Construct input chunks
        in_files = tqdm(list(in_dir.rglob('*')), desc='Collecting input files')
        inputs = sorted([
            path for path in tqdm(in_files, desc=f"Filtering input files with pattern \'{in_file_pattern}\'")
            if not in_file_pattern or re.match(in_file_pattern, path.name)
        ])
        if resume:
            out_files = list(self.out_dir.rglob('*.pdb')) \
                      + list(self.out_dir.rglob(f'*{NOPPI_EXTENSION}'))
            out_files = tqdm(out_files, desc='Collecting processed files')
            processed_pdbs = set(map(path_to_pdb_id, out_files))
            inputs = list(filter(lambda path: path_to_pdb_id(path) not in processed_pdbs, inputs))
        chunks = list(list_to_chunks(inputs, self.chunk_size))
        chunks = get_partition(chunks, *partition)

        # Run chunk processing in parallel
        if self._requires_bsa():
            self.dr_sasa = DR_SASA(
                tmp_dir=Path(f'./.dr_sasa_tmp_dir_{partition[0]}_{partition[1]}_{random_id(20)}'), 
                auto_clean='instant'  # Not to run out of disk memory
            )
        # with ProcessPoolExecutor(self.max_workers, mp_context=mp.get_context('fork')) as executor:
        with ProcessPoolExecutor(self.max_workers) as executor:
            future_to_chunk = {
                executor.submit(self._extract_chunk_paths, chunk): chunk 
                for chunk in chunks
            }
            futures = tqdm(as_completed(future_to_chunk), total=len(chunks))
            for future in futures:
                chunk = future_to_chunk[future]
                try:
                    future.result()
                except Exception as exc:
                    print(f'{chunk} generated an exception: {exc}')
                    print(traceback.format_exc(), end='\n\n')
        if self._requires_bsa():
            self.dr_sasa.clean()

    def _extract_chunk_paths(self, pdb_paths: Iterable[Path]) -> None:
        for pdb_path in pdb_paths:
            if self.verbose:
                print(f'[{datetime.datetime.now()}] {os.getpid()}: Starting {pdb_path}')
            self.extract(pdb_path)
            if self.verbose:
                print(f'[{datetime.datetime.now()}] {os.getpid()}: Finished {pdb_path}')

    def _preprocess_atom_df(
        self,
        atom_df: pd.DataFrame,
        chains: Optional[Iterable[str]] = None
    ) -> pd.DataFrame:
        """Preprocess atom dataframe in-place. Sorts residue identifiers and joins them into
        `res_ids`column. Applies filtering for protein chains. If `chains` is not only the 
        corresponding residues will be left.
        """
        # Rename columns from HADDOCK format
        if self.input_format == 'haddock':
            atom_df.chain_id = atom_df.segment_id
            atom_df.segment_id = ''

        # Filter chains
        if chains is not None:
            atom_df = atom_df[atom_df['chain_id'].isin(chains)]

        # Filter for 20 amino acid residues only.
        atom_df = atom_df[atom_df['residue_name'].isin(BASE_AMINO_ACIDS_3)]

        # Preapare residue ids
        node_id_columns = ['chain_id', 'residue_number', 'insertion']
        atom_df = atom_df.sort_values(node_id_columns)
        atom_df['res_id'] = atom_df.apply(
            lambda row: Residue(*[row[col] for col in node_id_columns]),
            axis=1
        )

        # Convert to float for memory efficiency
        for coord_col in ['x_coord', 'y_coord', 'z_coord']:
            atom_df[coord_col] = atom_df[coord_col].astype(np.float32)

        return atom_df

    def _get_radius_adjacency(
            self,
            atom_df: pd.DataFrame,
            radius: float
        ) -> pd.DataFrame:
        """Get residue-residue adjacency matrix based on atom contacts within predefined radius.
        """
        # Get adjacency matrix of atoms
        dists = sklearn.metrics.pairwise_distances(atom_df[['x_coord', 'y_coord', 'z_coord']])
        atom_adj = dists < radius

        # Get adjacency matrix of residues (at least one atom contact)
        res_ids = atom_df['res_id']
        atom_adj = pd.DataFrame(atom_adj, index=res_ids, columns=res_ids)
        res_adj = atom_adj.groupby(level=0).sum().T.groupby(level=0).sum() > 0
        return res_adj

    def _get_intra_adjacency(self, atom_df: pd.DataFrame) -> pd.DataFrame:
        """Get residue-residue adjacency matrix with True for all residues in the same chain.
        """
        # Construct adjacency matrix
        n_unique_rows = lambda df: len(df.drop_duplicates())
        n_res_per_chain = atom_df. \
            groupby('chain_id')[['residue_number', 'insertion']]. \
            apply(n_unique_rows)
        intra_adj = scipy.linalg.block_diag(*[np.ones((n, n), dtype=bool) for n in n_res_per_chain])

        # Wrap into annotated dataframe
        res_ids = atom_df['res_id'].unique()
        intra_adj = pd.DataFrame(intra_adj, index=res_ids, columns=res_ids)
        return intra_adj

    def _get_interacting_residues(
            self,
            adj: pd.DataFrame,
            a: str,
            b: str,
            symmetric_adj: bool = True
        ) -> set[Residue]:
        """Get interacting residues between chains `a` and `b` based on residue-residue adjacency
        matrix `adj`.
        """
        interacting = []
        for receptor, ligand in ((a, b), (b, a)) if not symmetric_adj else ((a, b),): 
            row_idx = pd.Series(adj.index).apply(lambda r: r[0] == receptor).to_numpy()
            col_idx = pd.Series(adj.columns).apply(lambda r: r[0] == ligand).to_numpy()
            sub_adj = adj.loc[row_idx, col_idx]
            side = sub_adj.replace(False, np.nan).stack().reset_index(allow_duplicates=True)
            side = side.iloc[:, :2].to_numpy().flatten()
            interacting.extend(side)
        return set(interacting)

    def _write_to_pdb(
            self,
            interface: Iterable[Residue],
            atom_df: pd.DataFrame,
            path: Union[os.PathLike, str],
            stats: dict = None,
            extra_columns: Sequence[str] = ('res_id', 'model_id')
        ) -> None:
        """Write atom data frame to PDB file.
        """
        # Filter atom_df for interfacial residues
        interface_atom_df = atom_df[atom_df['res_id'].isin(interface)]

        # Write to .pdb via Biopandas
        ppdb = PandasPdb()
        interface_atom_df = interface_atom_df.drop(columns=list(extra_columns))
        interface_atom_df['line_idx'] = np.arange(len(interface_atom_df))
        ppdb.df['ATOM'] = interface_atom_df
        self._add_extraction_remark(ppdb, stats=stats)
        ppdb.to_pdb(path)

    def _add_extraction_remark(self, ppdb: PandasPdb, code: int = 3, stats: dict = None) -> None:
        """Add REMARK entries on extracted PPI.
        """
        ppdb.add_remark(
            code,
            (f'PROTEIN-PROTEIN INTERACTION EXTRACTED WITH {PPIREF_NAME}'
            f' ({PPIREF_URL}).\n\n'
            f'EXTRACTION CRITERIA.')
        )
        ppdb.add_remark(code, f'KIND: {self.kind}', 1)
        ppdb.add_remark(code, f'EXTRACTION RADIUS: {self.radius} A', 1)
        ppdb.add_remark(code, f'EXPANSION RADIUS: {self.expansion_radius} A', 1)
        if stats is not None:
            ppdb.add_remark(code, '\nSTATISTICS.')
            for name, val in stats.items():
                ppdb.add_remark(code, f'{name}: {val}', 1)

    def _join_interfaces(
        self,
        interfaces: dict[tuple[str], set[Residue]]
    ) -> dict[tuple[str], set[Residue]]:
        raise NotImplementedError()
    
    def _requires_bsa(self) -> bool:
        return self.bsa or self.kind == 'bsa'
