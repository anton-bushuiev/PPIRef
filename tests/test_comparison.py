import shutil
from pathlib import Path

import numpy as np
import sklearn

from ppiref.comparison import IDist


def test_idist_deduplicate_embeddings():
    # Create redundant embeddings
    dummpy_tmp_dir = Path('./test_idist_deduplicate_embeddings')
    idist = IDist(pdb_dir=dummpy_tmp_dir, near_duplicate_threshold=0.04)
    idist.embeddings = {
        'A1': np.array([0.0, 0.0, 0.0]),
        'A2': np.array([0.0, 0.03, 0.0]),
        'B1': np.array([1.0, 0.0, 0.0]),
        'B2': np.array([1.0, 0.0, 0.0]),
        'C': np.array([2.0, 0.0, 0.0]),
    }
    for i in range(1000):
        idist.embeddings[f'D{i}'] = np.array([3.0, 0.0, 0.0])

    # Deduplicate
    sklearn.set_config(working_memory=0.5)
    idist.deduplicate_embeddings()

    # Clean
    shutil.rmtree(dummpy_tmp_dir)

    # Test no duplicates
    assert sorted(list(map(lambda x: x[0], idist.get_embeddings().index))) == ['A', 'B', 'C', 'D']
