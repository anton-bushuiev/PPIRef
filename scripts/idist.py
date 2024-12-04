import os
import warnings
import random
import pickle
from pathlib import Path

import click
import graphein

from ppiref.comparison import IDist
from ppiref.split import read_fold, read_split_source
from ppiref.definitions import PPIREF_DATA_DIR
from ppiref.utils.misc import get_partition

"""
Example run from the parent directory of this script:
python idist.py \
    --pairs="/scratch/project/open-29-57/antonb/idist2-private/data/benchmark/pairs_contrastive_6A_50k_pdbs.pkl" \
    --out_dir="/scratch/project/open-29-57/antonb/idist2-private/data/benchmark/idist" \
    --max_workers=126 \
    --partition_beg=0.00 \
    --partition_end=0.05
"""


graphein.verbose(enabled=False)
warnings.simplefilter(action='ignore', category=FutureWarning)


@click.command()
@click.option('--split', type=str, default='ppiref_6A_filtered')
@click.option('--fold', type=str, default='whole')
@click.option('--pairs', type=Path, default=None, help='Path to a pickled file with pairs of PPI paths to compare.')
@click.option('--pdb_dir', type=str, default=None)
@click.option('--out_dir', type=Path, default=None)
@click.option('--max_workers', type=int, default=os.cpu_count() - 2)
@click.option('--partition_beg', type=float, default=0.)
@click.option('--partition_end', type=float, default=1.)
@click.option('--seed', type=int, default=None)
@click.option('--out_file_suff', type=str, default='')
@click.option('--chunksize', type=int, default=1)
@click.option('--verbose', is_flag=True)
def main(split, fold, pairs, pdb_dir, out_dir, max_workers, partition_beg, partition_end, seed,
    out_file_suff, chunksize, verbose
):
    # Set seed (may be needed if `fold` is an integer number os random samples)
    if seed is not None:
        random.seed(seed)

    # Read paths to all PPIs and partition them
    if pairs is not None:  # Pairs provided
        name = f'{pairs.stem}'
        with open(pairs, 'rb') as f:
            ppis0, ppis1 = pickle.load(f)
        pairs = list(zip(ppis0, ppis1))
        pairs = get_partition(pairs, partition_beg, partition_end)
        inputs = {'ppi_pairs': pairs}
        if verbose:
            print(f'Read {len(pairs)} PPI pairs')
    else:  # Split and fold provided
        name = f'{split},{fold}'
        ppis0 = read_fold(split, fold)
        ppis1 = read_fold(split, fold)
        ppis0 = get_partition(ppis0, partition_beg, partition_end)
        ppis1 = get_partition(ppis1, partition_beg, partition_end)
        inputs = {'ppis0': ppis0, 'ppis1': ppis1}
        if verbose:
            print(f'Read {len(ppis0)} and {len(ppis1)} PPIs')

    # Init path to directory with full PDBs
    if pdb_dir is None:
        pdb_dir = PPIREF_DATA_DIR / 'pdb'

    # Init path to results
    if out_dir is None:
        out_dir = read_split_source(split).parent / 'clustering'
    out_dir.mkdir(exist_ok=True, parents=True)
    print(f'Results will be saved to {out_dir}')

    # Init comparator
    idist = IDist(pdb_dir=pdb_dir, max_workers=max_workers, verbose=verbose)

    # Embed
    ppis = list(set(ppis0 + ppis1)) if pairs is not None else ppis0
    idist.embed_parallel(ppis, chunksize=chunksize)
    df_embeddings = idist.get_embeddings()
    df_embeddings.to_csv(out_dir / f'idist_emb_{name}_{partition_beg}-{partition_end}{out_file_suff}.csv')

    # Compare
    df_comp = idist.compare_all_against_all(**inputs, embed=False)
    df_comp.to_csv(out_dir / f'idist_distmat_{name}_{partition_beg}-{partition_end}{out_file_suff}.csv', index=False)


if __name__ == '__main__':
    main()
