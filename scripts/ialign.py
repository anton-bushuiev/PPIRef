import os
import random
from pathlib import Path

import click

from ppiref.comparison import IAlign
from ppiref.split import read_fold, read_split_source


# Example:
# python ialign.py --split="skempi2_cleaned_rde_net" --fold="whole" --max_workers=64

@click.command()
@click.option('--split', type=str, default='ppiref_dips_filtered_sample_500')
@click.option('--fold', type=str, default='whole')
@click.option('--out_dir', type=str, default=None)
@click.option('--max_workers', type=int, default=os.cpu_count() - 2)
@click.option('--seed', type=int, default=None)
@click.option('--out_file_suff', type=str, default='')
@click.option('--verbose', '-v', is_flag=True)
def main(split, fold, out_dir, max_workers, seed, out_file_suff, verbose):
    # Set seed (may be needed if `fold` is an integer number of random samples)
    if seed is not None:
        random.seed(seed)

    # Init paths to all PPIs
    ppis0 = read_fold("skempi2_complexes_covid_sak", "whole")
    ppis1 = read_fold("skempi2_complexes_covid_sak", "whole")

    # Init path to results directory
    if out_dir is None:
        out_dir = read_split_source(split).parent / 'clustering'
    out_dir.mkdir(exist_ok=True, parents=True)

    # Init comparator
    ialign = IAlign(max_workers=max_workers, verbose=verbose)

    # Compare
    df = ialign.compare_all_against_all(ppis0, ppis1)
    out_file = out_dir / f'ialign_comp_{split},{fold}{out_file_suff}.csv'
    if verbose:
        print(f'Result data frame ({len(df)} lines) saved to {out_file}')
    df.to_csv(out_file)


if __name__ == '__main__':
    main()
