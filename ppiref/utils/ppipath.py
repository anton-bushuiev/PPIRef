"""PPIPath class and related functions for PPI file paths."""
from typing import Union, Iterable
from pathlib import Path


NOPPI_EXTENSION: str = '.noppi'


# TODO Check for the absence of '_' symbols in input paths
class PPIPath(type(Path())):
    """Path to a protein-protein interaction (PPI) file."""

    def __new__(cls, *args, **kwargs):
        """https://stackoverflow.com/q/61689391"""
        obj = super().__new__(cls, *args, **kwargs)
        if obj.suffix not in ['.pdb', NOPPI_EXTENSION]:
            raise ValueError(f'Wrong PPI file extension {obj.suffix}.')
        return obj        

    def is_dummy(self) -> bool:
        """Check if the PPI file is a dummy file that does not contain any PPIs."""
        return self.suffix == NOPPI_EXTENSION
    
    def ppi_id(self) -> str:
        """Get PPI ID from the path (1abc_A_B)."""
        return path_to_ppi_id(self)

    def pdb_id(self) -> str:
        """Get PDB ID from the path (1abc)."""
        return path_to_pdb_id(self)

    def partners(self) -> list[str]:
        """Get the list of partners from the path (['A', 'B'])."""
        return path_to_partners(self)
    
    def is_equivalent_to(self, other: Union[Path, str]) -> bool:
        """Check if two PPI paths are equivalent, up to partners order."""
        other = PPIPath(other)
        return (
            self.pdb_id() == other.pdb_id() and
            set(self.partners()) == set(other.partners())
        )

    def contains(self, other: Union[Path, str]) -> bool:
        """Check if one PPI path contains another based on sets of partners."""
        other = PPIPath(other)
        return (
            self.pdb_id() == other.pdb_id() and
            set(self.partners()) >= set(other.partners())
        )

    def with_sorted_partners(self) -> Path:
        """Get a path with partners sorted alphabetically."""
        return self.with_stem(self.pdb_id() + '_' + '_'.join(sorted(self.partners())))
    

    @classmethod
    def construct(cls, pref: Union[str, Path], name: str, partners: str, nested: bool = True) -> 'PPIPath':
        """Construct a PPI path from the components."""
        pref = Path(pref)
        partners = '_'.join(partners)
        if nested:
            pref /= name[1:3]
        return cls(pref / f'{name}_{partners}.pdb')

def path_to_ppi_id(path: Union[Path, str]) -> str:
    """Get PPI ID from the path (1abc_A_B)."""
    stem = Path(path).stem
    if stem[0] == '.':  # hidden file
        stem = stem[1:]
    return stem.split('.', 1)[0]


def path_to_pdb_id(path: Union[Path, str]) -> str:
    """Get PDB ID from the path (1abc)."""
    return path_to_ppi_id(path).split('_')[0]


def path_to_partners(path: Union[Path, str]) -> list:
    """Get the list of partners from the path (['A', 'B'])."""
    return path_to_ppi_id(path).split('_')[1:]


def pdb_id_to_nested_suffix(pdb_id: str) -> str:
    """Get a nested suffix for a PDB ID (1abc -> ab/1abc.pdb)."""
    return f'{pdb_id[1:3]}/{pdb_id}.pdb'


def ppi_id_to_nested_suffix(ppi_id: str) -> str:
    """Get a nested suffix for a PPI ID (1abc_A_B -> ab/1abc_A_B.pdb)."""
    return pdb_id_to_nested_suffix(ppi_id)
