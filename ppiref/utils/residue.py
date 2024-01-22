from collections import namedtuple

from Bio.Data.PDBData import protein_letters_1to3
from graphein.protein.resi_atoms import BASE_AMINO_ACIDS


# Amino acid codes
BASE_AMINO_ACIDS = BASE_AMINO_ACIDS
BASE_AMINO_ACIDS_INVERSE = {
    aa: i for i, aa in enumerate(BASE_AMINO_ACIDS)
}
BASE_AMINO_ACIDS_GROUPED = [
    'R', 'H', 'K',
    'D', 'E',
    'S', 'T', 'N', 'Q',
    'C', 'G', 'P',
    'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'
]
BASE_AMINO_ACIDS_3 = list(map(
    lambda x: protein_letters_1to3[x], BASE_AMINO_ACIDS
))


# Helper namedtuple to store unique residue ids
Residue = namedtuple(
    'Residue',
    ['chain_id', 'residue_number', 'insertion'],
    defaults=(None, None, '')
)
