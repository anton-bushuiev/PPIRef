"""Example
python scripts/extract_ppis_in_parallel.py \
    --in_dir=/scratch/project/open-26-23/antonb/PPIRef/data/pdb \
    --out_dir=/scratch/project/open-26-23/antonb/PPIRef/data/ppiref/ppi
"""
import os
from pathlib import Path

import click

from ppiref.extraction import PPIExtractor
from ppiref.definitions import *


@click.command()
@click.option('--in_dir', type=Path)
@click.option('--out_dir', type=Path)
@click.option('--kind', type=str, default='heavy')
@click.option('--radius', type=float, default=10.)
@click.option('--expansion_radius', type=float, default=0.)
@click.option('--bsa', type=bool, default=True)
@click.option('--join', type=bool, default=False)
@click.option('--nest_out_dir', type=bool, default=True)
@click.option('--max_workers', type=int, default=os.cpu_count() - 2)
@click.option('--chunk_size', type=int, default=1)
@click.option('--partition_beg', type=float, default=0.)
@click.option('--partition_end', type=float, default=1.)
@click.option('--input_format', type=str, default="pdb")
def main(
    in_dir, out_dir, kind, radius, expansion_radius, bsa, join,
    nest_out_dir, max_workers, chunk_size, partition_beg, partition_end, input_format
):
    ppi_extractor = PPIExtractor(
        out_dir=out_dir,
        kind=kind,
        radius=radius,
        expansion_radius=expansion_radius,
        bsa=bsa,
        join=join,
        nest_out_dir=nest_out_dir,
        max_workers=max_workers,
        chunk_size=chunk_size,
        verbose=True,
        input_format=input_format
    )
    ppi_extractor.extract_parallel(
        in_dir=in_dir,
        partition=(partition_beg, partition_end)
    )


if __name__ == '__main__':
    main()
