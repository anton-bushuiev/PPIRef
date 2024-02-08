import os
import warnings
import random

import click
import graphein

from ppiref.comparison import IDist
from ppiref.split import read_fold, read_split_source
from ppiref.definitions import PPIREF_DATA_DIR


graphein.verbose(enabled=False)
warnings.simplefilter(action='ignore', category=FutureWarning)


@click.command()
@click.option('--split', type=str, default='ppiref_6A_filtered')
@click.option('--fold', type=str, default='whole')
@click.option('--pdb_dir', type=str, default=None)
@click.option('--out_dir', type=str, default=None)
@click.option('--max_workers', type=int, default=os.cpu_count() - 2)
@click.option('--seed', type=int, default=None)
@click.option('--out_file_suff', type=str, default='')
def main(split, fold, pdb_dir, out_dir, max_workers, seed, out_file_suff):
    # Set seed (may be needed if `fold` is an integer number os random samples)
    if seed is not None:
        random.seed(seed)

    # Init paths to all PPIs
    ppis = read_fold(split, fold)

    # Init path to directory with full PDBs
    if pdb_dir is None:
        pdb_dir = PPIREF_DATA_DIR / 'pdb'

    # Init path to results directory
    if out_dir is None:
        out_dir = read_split_source(split).parent / 'clustering'
    out_dir.mkdir(exist_ok=True, parents=True)

    # Init comparator
    idist = IDist(pdb_dir=pdb_dir, max_workers=max_workers, verbose=True)

    # Embed
    idist.embed_parallel(ppis)
    df_embeddings = idist.get_embeddings()
    df_embeddings.to_csv(out_dir / f'idist_emb_{split},{fold}{out_file_suff}.csv')

    # Compare
    df_comp = idist.compare_all_against_all(ppis, ppis, embed=False)
    df_comp.to_csv(out_dir / f'idist_distmat_{split},{fold}.csv', index=False)


if __name__ == '__main__':
    main()
