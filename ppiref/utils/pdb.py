import warnings
import urllib
from typing import Optional, Union
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser
from biopandas.pdb import PandasPdb


def get_n_models(pdb_path: Path) -> int:
    pdb = PDBParser(QUIET=True)
    structure = pdb.get_structure('', pdb_path)
    return len(list(structure))


def get_first_model(pdb_path: Path, out_path: Optional[Path] = None) -> None:
    # Parse to first model in biopandas
    ppdb_df = PandasPdb().read_pdb(str(pdb_path))
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', category=pd.errors.SettingWithCopyWarning)
        ppdb_df = ppdb_df.get_model(1)

    # Write to .pdb
    if out_path is not None:
        del ppdb_df.df['ATOM']['model_id']
        ppdb_df.to_pdb(str(out_path), records=['ATOM'])  # TODO other records?

    return ppdb_df


def download_pdb(
    pdb_id: str,
    dir: Optional[Union[Path, str]] = '.',
    path: Optional[Path] = None
) -> None:
    if path is None:
        path = Path(dir) / f'{pdb_id}.pdb'
    urllib.request.urlretrieve(
        f'http://files.rcsb.org/download/{pdb_id}.pdb',
        str(path)
    )
