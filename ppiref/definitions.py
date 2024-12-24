"""Global variables used across the package."""
import pathlib

PPIREF_NAME = 'PPIRef'
PPIREF_URL = 'https://github.com/anton-bushuiev/PPIRef'

# Dirs
PPIREF_ROOT_DIR = pathlib.Path(__file__).parent.absolute()
PPIREF_REPO_DIR = PPIREF_ROOT_DIR.parent
PPIREF_DATA_DIR = PPIREF_ROOT_DIR / 'data'
PPIREF_SPLITS_DIR = PPIREF_DATA_DIR / 'splits'
PPIREF_TEST_DATA_DIR = PPIREF_DATA_DIR / 'test'
PPIREF_EXTERNAL_DIR = PPIREF_REPO_DIR / 'external'

# BM5
PPIREF_BM5_CSV_PATH = PPIREF_DATA_DIR / 'bm5/bm5_merge.csv'

# External software
DR_SASA_PATH = PPIREF_EXTERNAL_DIR / 'dr_sasa_n/build/dr_sasa' 
IALIGN_PATH = PPIREF_EXTERNAL_DIR / 'ialign/bin/ialign.pl'
USALIGN_PATH = PPIREF_EXTERNAL_DIR / 'USalign/USalign'
FOLDSEEK_PATH = PPIREF_EXTERNAL_DIR / 'foldseek/bin/foldseek'
