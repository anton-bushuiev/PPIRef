from ppiref.surface import DR_SASA
from ppiref.definitions import PPIREF_TEST_DATA_DIR


def test_dr_sasa_parse_residue():
    parsed = DR_SASA.parse_residue('VAL/M/14A')
    assert tuple(parsed) == ('M', 14, 'A')

    parsed = DR_SASA.parse_residue('VAL/M/141')
    assert tuple(parsed) == ('M', 141, '')

    parsed = DR_SASA.parse_residue('VAL/M/-141')
    assert tuple(parsed) == ('M', -141, '')


def test_dr_sasa_nmr():
    pdb_path = PPIREF_TEST_DATA_DIR / '1a0n.pdb'
    dr_sasa = DR_SASA()
    buried_residues, bsa = dr_sasa(pdb_path, ('A', 'B'))
    assert len(set(buried_residues)) == len(buried_residues)
    assert bsa < 500


# TODO Test BSA from full file is same as from extracted PPI file
