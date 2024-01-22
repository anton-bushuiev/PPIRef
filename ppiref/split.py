import json
import warnings
import itertools
import random
from tqdm import tqdm

from pathlib import Path
from typing import Iterable, Union

from ppiref.utils.ppipath import PPIPath, ppi_id_to_nested_suffix
from ppiref.definitions import PPIREF_ROOT_DIR, PPIREF_SPLITS_DIR


# TODO `Split` class (probably `Fold` as well)


def write_split(
    location: Union[Path, str],
    source: Path,
    folds: dict[str, Iterable[Union[str, Path]]],
    analyze: bool = True
) -> None:
    # Process args
    location = _process_location(location)
    for fold in folds:
        folds[fold] = [
            Path(p).stem if (isinstance(p, Path) or '/' in p) else p
            for p in folds[fold]
        ]
    
    # Analyze folds
    if analyze:
        # PPIs exist
        for fold, ppis in folds.items():
            for ppi in ppis:
                path = source / ppi_id_to_nested_suffix(ppi)
                if not path.is_file():
                    warnings.warn(f'{path} from fold \'{fold}\' does not exist.')

        # names unique withing folds
        for fold, ppis in folds.items():
            if len(ppis) != len(set(ppis)):
                warnings.warn(f'Repetitive PPI names in \'{fold}\'.')

        # all folds disjoint
        fold_pairs = itertools.combinations(folds.items(), 2)
        for (fold_a, ppis_a), (fold_b, ppis_b) in fold_pairs:
            if set(ppis_a) & set(ppis_b):
                warnings.warn(f'Folds {fold_a} and {fold_b} are not disjoint.')

        # cover complete PPI source
        ppis_folds = set().union(*folds.values())
        ppis_source = set([path.stem for path in source.rglob('*.pdb')])
        if len(ppis_folds) != len(ppis_source):
            warnings.warn(
                f'Split is not complete: {len(ppis_folds)} of {len(ppis_source)} PPIs contained.'
        )

    # Write to .json
    source = source.relative_to(PPIREF_ROOT_DIR)
    split = dict(source=str(source), folds=folds)
    json_object = json.dumps(split)
    with open(location, 'w') as file:
        file.write(json_object)


def read_split(
    location: Union[Path, str],
    full_paths: bool = True
) -> dict[str, list[Path]]:
    location = _process_location(location)

    # Read .json
    with open(location, 'r') as file:
        split = json.load(file)

    # Parse source
    source = Path(split['source'])
    if not source.is_absolute():
        source = PPIREF_ROOT_DIR / source

    # Parse folds
    folds = split['folds']
    if full_paths:
        for fold in tqdm(folds, desc='converting paths format', leave=False):
            folds[fold] = [PPIPath(source / ppi_id_to_nested_suffix(p)) for p in folds[fold]]
    return folds
    

def read_fold(
    location: Union[Path, str],
    fold: Union[str, int],
    full_paths: bool = True,
    processed_split: dict[str, list[Path]] = None
) -> list[Path]:
    if processed_split is None:
        split = read_split(location, full_paths=full_paths)
    else:
        split = processed_split
    
    if fold in split:
        return split[fold]
    elif '+' in fold:
        return sum([read_fold(location, f, processed_split=split) for f in fold.split('+')], start=[])
    elif fold == 'whole':
        return sum(split.values(), start=[])
    elif (isinstance(fold, int) and fold >= 0) or fold.isdigit():  # non-negative int
        return random.sample(read_fold(location, 'whole'), int(fold))
    else:
        raise ValueError(f'Wrong `fold` value {fold}.')
    

def read_split_source(
    location: Union[Path, str]
) -> Path:
    location = _process_location(location)

    # Read .json
    with open(location, 'r') as file:
        split = json.load(file)

    # Parse source
    source = Path(split['source'])
    if not source.is_absolute():
        source = PPIREF_ROOT_DIR / source

    return source


def _process_location(location: Union[Path, str]) -> Path:
    if not isinstance(location, str) and not isinstance(location, Path):
        raise ValueError('Wrong `location` type.')
    elif isinstance(location, str):
        if location.startswith('/'):
            location = Path(location)
        else:
            location = (PPIREF_SPLITS_DIR / location).with_suffix('.json')
    else:  # isinstance(location, Path) 
        if location.suffix != '.json':
            raise ValueError(f'Wrong `location` extension {location.suffix}.')
    return location
