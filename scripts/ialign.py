import os
import random
import pickle
from pathlib import Path

import click

from ppiref.comparison import IAlign
from ppiref.split import read_fold, read_split_source
from ppiref.utils.misc import get_partition

# Examples:
# python ialign.py --split="skempi2_cleaned_rde_net" --fold="whole" --max_workers=64

@click.command()
@click.option('--split', type=str, default='ppiref_dips_filtered_sample_500')
@click.option('--fold', type=str, default='whole')
@click.option('--pairs', type=Path, default=None, help='Path to a pickled file with pairs of PPI paths to compare.')
@click.option('--out_dir', type=Path, default=None)
@click.option('--max_workers', type=int, default=os.cpu_count() - 2)
@click.option('--partition_beg', type=float, default=0.)
@click.option('--partition_end', type=float, default=1.)
@click.option('--seed', type=int, default=None)
@click.option('--out_file_suff', type=str, default='')
@click.option('--verbose', '-v', is_flag=True)
def main(split, fold, pairs, out_dir, max_workers, partition_beg, partition_end, seed,
    out_file_suff, verbose
):
    # Set seed (may be needed if `fold` is an integer number of random samples)
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

    # Init path to results directory
    if out_dir is None:
        out_dir = read_split_source(split).parent / 'clustering'
    out_dir.mkdir(exist_ok=True, parents=True)
    if verbose:
        print(f'Results will be saved to {out_dir}')

    # Init comparator
    ialign = IAlign(max_workers=max_workers, verbose=verbose)
    if verbose:
        print(f'Initialized comparator with {max_workers} workers')

    # Compare
    df = ialign.compare_all_against_all(**inputs)
    out_file = out_dir / f'ialign_comp_{name}_{partition_beg}-{partition_end}{out_file_suff}.csv'
    if verbose:
        print(f'Result data frame ({len(df)} lines) saved to {out_file}')
    df.to_csv(out_file, index=False)


if __name__ == '__main__':
    main()
