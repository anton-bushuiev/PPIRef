import shutil

import dill
import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb

from ppiref.definitions import PPIREF_TEST_DATA_DIR
from ppiref.extraction import PPIExtractor
from ppiref.utils.residue import Residue
from ppiref.utils.ppipath import path_to_partners


def test_dips_reproduction():
    """Test that extraction matches to the DIPS extraction done with atom3.
    Note: in DIPS data frames insertion codes are merged with the string `residue` column while
    in PPIRef the raw PDBs are stored and residue numbers and insertion codes are read separately
    into Biopandas dataframe. This function does not test insertion codes.
    """
    path_pdb = PPIREF_TEST_DATA_DIR / '10gs.pdb'
    path_dips_pair = PPIREF_TEST_DATA_DIR / '10gs.pdb1_0.dill'

    # Extract interface
    out_dir = PPIREF_TEST_DATA_DIR / 'tmp'
    ppi_extractor = PPIExtractor(kind='heavy', radius=6., out_dir=out_dir)
    ppi_extractor.extract(path_pdb)
    df_ppiref = PandasPdb().read_pdb(str(PPIREF_TEST_DATA_DIR / 'tmp/0g/10gs_A_B.pdb')).df['ATOM']
    interface = df_ppiref[['chain_id', 'residue_number']].to_numpy().tolist()
    interface = set(map(lambda x: Residue(*x, ''), interface))

    # Read DIPS interface pair
    with open(PPIREF_TEST_DATA_DIR / '10gs.pdb1_0.dill', 'rb') as file:
        dips_pair = dill.load(file)
    dfs = [dips_pair.df0, dips_pair.df1]
    for i in range(2):
        pos_res = dfs[i].loc[np.unique(dips_pair.pos_idx[:, i])]['residue']
        dfs[i] = dfs[i][dfs[i]['residue'].isin(pos_res)]
    df_dips = pd.concat(dfs)
    df_dips['residue'] = df_dips['residue'].astype(int)
    interface_dips = df_dips[['chain', 'residue']].to_numpy().tolist()
    interface_dips = set(map(lambda x: Residue(*x, ''), interface_dips))

    # Compare to DIPS
    assert interface == interface_dips

    # Clean
    shutil.rmtree(out_dir)


def test_no_nucleic_acids():
    """Test that nucleic acids are not considered for interfaces.
    """
    # Extact PPIs from complex of DNA and protein mixture
    out_dir = PPIREF_TEST_DATA_DIR / 'tmp'
    ppi_extractor = PPIExtractor(out_dir=out_dir, radius=10.)
    ppi_extractor.extract(PPIREF_TEST_DATA_DIR / '1a02.pdb')

    # Test no nucleic acid chains in extracted interfaces
    ppi_paths = list((PPIREF_TEST_DATA_DIR / 'tmp').rglob('*.pdb'))
    assert len(ppi_paths)
    for path in ppi_paths:
        partners = path_to_partners(path)
        assert 'A' not in partners and 'B' not in partners

    # Clean
    shutil.rmtree(out_dir)


def test_bsa_no_transitive_interactions():
    """Test that the BSA-based interface extraction does not suffer from "transitive interactions"
    unlike radius-based extraction. The tested pdb contains 7 chains with a cyclic symmetry of
    contacts and therefore should contain 7 interactions. Heavy 10A extraction generates 21 PPIs.
    """
    pdb_path = PPIREF_TEST_DATA_DIR / '1aon_OPQRSTU.pdb'
    out_dir = PPIREF_TEST_DATA_DIR / 'tmp_test_bsa_no_transitive_interactions'
    ppi_extractor = PPIExtractor(out_dir, kind='bsa', radius=6.)
    ppi_extractor.extract(pdb_path)
    n_ppis = len(list(out_dir.rglob('*.pdb')))
    shutil.rmtree(out_dir)
    assert n_ppis  == 7


def test_bsa_is_radius_subset():
    """Test that the interface defined by buried residues is a subset of the radius-based one with
    sufficiently big raidus.
    """
    # Extract
    pdb_path = PPIREF_TEST_DATA_DIR / '1a02.pdb'
    out_dir_bsa = PPIREF_TEST_DATA_DIR / 'tmp_bsa'
    out_dir_radius = PPIREF_TEST_DATA_DIR / 'tmp_radius'
    ppi_extractor_bsa = PPIExtractor(out_dir_bsa, kind='bsa', radius=10.)
    ppi_extractor_radius = PPIExtractor(out_dir_radius, kind='heavy', radius=10.)
    ppi_extractor_bsa.extract(pdb_path)
    ppi_extractor_radius.extract(pdb_path)

    # Read and test for subsets
    for (a, b) in (('F', 'J'), ('F', 'N'), ('J', 'N')):
        ppi_path_suffix = f'a0/1a02_{a}_{b}.pdb'
        df_bsa = PandasPdb().read_pdb(str(out_dir_bsa / ppi_path_suffix)).df['ATOM']
        df_radius = PandasPdb().read_pdb(str(out_dir_radius / ppi_path_suffix)).df['ATOM']
        assert len(df_bsa) < len(df_radius)
        assert set(df_bsa['atom_number']) < set(df_radius['atom_number'])

    # Clean
    shutil.rmtree(out_dir_bsa)
    shutil.rmtree(out_dir_radius)


def test_expansion_hierarchy():
    """Test that interfaces with gradually increasing expansion radiuses are inclusively ordered.
    """
    pdb_path = PPIREF_TEST_DATA_DIR / '10gs.pdb'
    atoms = []
    for r in (0., 4., 10.):
        out_dir = PPIREF_TEST_DATA_DIR / 'tmp_test_expansion_hierarchy'
        ppi_extractor = PPIExtractor(out_dir, kind='heavy', radius=4., expansion_radius=r)
        ppi_extractor.extract(pdb_path)
        df = PandasPdb().read_pdb(str(out_dir / '0g/10gs_A_B.pdb')).df['ATOM']
        atoms.append(set(df['atom_number']))
    shutil.rmtree(out_dir)
    assert atoms[0] < atoms[1]
    assert atoms[0] < atoms[2]
    assert atoms[1] < atoms[2]


def test_partner_specification():
    """Test correct extraction of a PPI with specified chains.
    """
    pdb_path = PPIREF_TEST_DATA_DIR / '1ahw.pdb'
    partners = ['A', 'B', 'C']
    out_dir = PPIREF_TEST_DATA_DIR / 'tmp_test_partner_specification'

    # No joining
    ppi_extractor = PPIExtractor(out_dir)
    ppi_extractor.extract(pdb_path, partners=partners)
    n_ppis = len(list(out_dir.rglob('*.pdb')))
    shutil.rmtree(out_dir)
    assert n_ppis == 3    

    # Join
    ppi_extractor = PPIExtractor(out_dir, join=True)
    ppi_extractor.extract(pdb_path, partners=partners)
    ppi_paths = list(out_dir.rglob('*.pdb'))
    n_ppis = len(ppi_paths)
    suffix_partners = path_to_partners(ppi_paths[0])
    pdb_partners = PandasPdb().read_pdb(ppi_paths[0]).df['ATOM']['chain_id'].unique().tolist()
    shutil.rmtree(out_dir)
    assert n_ppis == 1
    assert suffix_partners == partners
    assert pdb_partners == partners    


def test_haddock_format():
    """Test that for a HADDOCK file the subset of atoms is extracted.
    """
    pdb_path = PPIREF_TEST_DATA_DIR / '1MAH-ti5-it0-805.pdb'
    out_dir = PPIREF_TEST_DATA_DIR / 'tmp_test_haddock_format'
    ppi_extractor = PPIExtractor(out_dir, radius=10., input_format='haddock')
    ppi_extractor.extract(pdb_path)

    df_interface = PandasPdb().read_pdb(str(out_dir / 'MA/1MAH-ti5-it0-805_A_B.pdb')).df['ATOM']
    df_complex = PandasPdb().read_pdb(str(pdb_path)).df['ATOM']
    
    shutil.rmtree(out_dir)
    assert set(df_interface['atom_number']) < set(df_complex['atom_number'])
