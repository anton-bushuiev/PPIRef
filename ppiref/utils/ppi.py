"""Utilities for working with and analyzing PPI files."""
from pathlib import Path
from typing import Any, Union

from biopandas.pdb import PandasPdb

from ppiref.visualization import PyMOL
from ppiref.utils.ppipath import PPIPath, ppi_id_to_nested_suffix
from ppiref.definitions import PPIREF_DATA_DIR


PYMOL = PyMOL()


class PPI:
    def __init__(self, ppi: Union[Path, str]) -> None:
        """This class represents a protein-protein interaction (PPI) defined by its path
        and provides methods for visualizing and analyzing it.

        Args:
            ppi (Union[Path, str]): Either full path or PPI id to get from PPIRef.
        """
        # Process args
        if isinstance (ppi, str) and '.' not in ppi:
            path = PPIREF_DATA_DIR / 'ppiref/ppi_6A' / ppi_id_to_nested_suffix(ppi)
        else:
            path = ppi
        self.path = PPIPath(path)
        if self.is_dummy():
            return

        # Parse atoms
        self.ppdb_df = PandasPdb().read_pdb(str(path))

        # Parse remarks into statistics
        remarks = self.ppdb_df.df['OTHERS']
        remarks = remarks[remarks['record_name'] == 'REMARK']['entry']
        remarks = remarks[remarks.str.contains(': ')].apply(self._parse_stat_remark)
        self.stats = dict(remarks.tolist())

    def is_dummy(self) -> bool:
        """Check if the PPI is a dummy file that does not contain any PPIs."""
        return self.path.is_dummy()

    def visualize(self, **kwargs):
        """Visualize the PPI in PyMOL."""
        return PYMOL.display_ppi(self.path, **kwargs)
    
    def n_residues(self) -> int:
        """Get the number of residues in the PPI."""
        return len(self.ppdb_df.df['ATOM'][['chain_id', 'residue_number', 'insertion']].drop_duplicates())
    
    def _parse_stat_remark(self, remark: str) -> tuple[str, Any]:
        """Parse a remark entry into a name-value pair."""
        name, val = remark.split(': ', maxsplit=1)

        # Name
        name = name.split(maxsplit=1)[1]

        # Value
        val = val.strip()
        if val.endswith('A') or val.endswith('A^2'):
            val = val.split()[0]
            if val == 'None':
                val = None
            else:
                val = float(val)

        return name, val


def sort_partners(ppi_id):
    """Sort partners in the ``ppi_id``."""
    parts = ppi_id.split('_')
    pdb_id, partners = parts[0], parts[1:]
    return pdb_id + '_' + '_'.join(sorted(partners))
