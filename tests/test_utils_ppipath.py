from ppiref.utils.ppipath import PPIPath


def is_equivalent_to():
    assert PPIPath('root/ab/1abc_A_E_D.pdb').is_equivalent_to('root/ab/1abc_A_D_E.pdb')
    assert not PPIPath('root/ab/1abc_A_E_D.pdb').is_equivalent_to('root/ab/1abc_A_D.pdb')
    assert PPIPath('root/ab/1abc_A_E_D.pdb').is_equivalent_to('root/ab/2abc_A_D_E.pdb')


def test_contains():
    assert PPIPath('root/ab/1abc_A_E_D.pdb').contains('root/ab/1abc_A_D.pdb')
    assert PPIPath('root/ab/1abc_A_E_D.pdb').contains('root/ab/1abc_A_E_D.pdb')
    assert PPIPath('root/ab/1abc_A_E_D.pdb').contains('root/ab/1abc_A.pdb')
    assert PPIPath('root/ab/1abc_A_E.pdb').contains('root/ab/1abc_A_E.pdb')
    assert not PPIPath('root/ab/1abc_A_E_D.pdb').contains('root/ab/1abc_A_K.pdb')
    assert not PPIPath('root/ab/1abc_A_E.pdb').contains('root/ab/1abc_A_E_D.pdb')
    assert not PPIPath('root/ab/1abc_A_E_D.pdb').contains('root/ab/2abc_A_D.pdb')


def test_with_sorted_partners():
    assert PPIPath('root/ab/1abc_A_E_D.pdb').with_sorted_partners() == PPIPath('root/ab/1abc_A_D_E.pdb')
