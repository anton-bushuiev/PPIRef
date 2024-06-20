"""Functions to read and write data splits in a standardized way and single JSON format."""
import json
import warnings
import itertools
import random
from tqdm import tqdm

from pathlib import Path
from typing import Iterable, Union

from ppiref.utils.ppipath import PPIPath, ppi_id_to_nested_suffix
from ppiref.definitions import PPIREF_ROOT_DIR, PPIREF_SPLITS_DIR


def write_split(
    location: Union[Path, str],
    source: Path,
    folds: dict[str, Iterable[Union[str, Path]]],
    analyze: bool = True
) -> None:
    """
    Write data split of protein-protein interactions to a JSON file.

    Args:
        location (Union[Path, str]): Destination of the split file. This can be:
        
            - A name of the split file, which will be stored as ``ppiref.definitions.PPIREF_SPLITS_DIR / f'{location}.json'``.
            - A Path object representing the path to the JSON split file.
            - An absolute path string starting with ``'/'`` representing the path to the JSON split file.
            
        source (Path): Path to the source directory containing PPIs.
        folds (dict[str, Iterable[Union[str, Path]]]): Dictionary of folds. Each fold (e.g. ``'train'``
            or ``'val'``) is a list of PPIs represented by their names or paths.
        analyze (bool, optional): If True, run simple sanity checks such as no overlapping PPI ids
            across the folds. If any issues are found, warnings are raised. Defaults to True.
    """
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
) -> dict[str, list[Union[Path, str]]]:
    """
    Read data split of protein-protein interactions from a JSON file.

    Args:
        location (Union[Path, str]): Source of the split file. This can be:
        
            - A name of the split file, which is stored as ``ppiref.definitions.PPIREF_SPLITS_DIR / f'{location}.json'``.
            - A Path object representing the path to the JSON split file.
            - An absolute path string starting with ``'/'`` representing the path to the JSON split file.
            
        full_paths (bool, optional): If set to True, return full paths instead of IDs. Defaults to True.

    Returns:
        dict[str, list[Union[Path, str]]]: Dictionary of data folds. Each fold is a list of PPIs 
        represented by their IDs or paths.
    """
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
) -> list[Union[Path, str]]:
    """Read a specific data fold from a data split of protein-protein interactions.

    Args:
        location (Union[Path, str]): Source of the split file. This can be:
        
            - A name of the split file, which is stored as ``ppiref.definitions.PPIREF_SPLITS_DIR / f'{location}.json'``.
            - A Path object representing the path to the JSON split file.
            - An absolute path string starting with ``'/'`` representing the path to the JSON split file.
            
        fold (Union[str, int]): Name of the fold to read. Should match one of the keys in the split
            dictionary. If ``'whole'``, all PPIs are returned. If a ``'+'``-separated string, PPIs from all
            specified folds are returned. If a non-negative int, a random sample of that size is
            returned. Defaults to ``'whole'``.
            
        full_paths (bool, optional): If set to True, return full paths instead of IDs. Defaults to True.
        
        processed_split (dict[str, list[Path]], optional): Pre-defined split dictionary. Defaults to None.

    Returns:
        list[Union[Path, str]]: Fold of PPIs represented by the list of their IDs or paths.
    """
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
    """
    Read source directory containing .pdb files in a data split of protein-protein interactions.

    Args:
        location (Union[Path, str]): Source of the split file. This can be:
        
            - A name of the split file, which is stored as ``ppiref.definitions.PPIREF_SPLITS_DIR / f'{location}.json'``.
            - A Path object representing the path to the JSON split file.
            - An absolute path string starting with ``'/'`` representing the path to the JSON split file.

    Returns:
        Path: Path to the source directory containing .pdb files with PPI structures.
    """
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
    """
    Convert data split (destination or source) location to a Path object.

    Args:
        location (Union[Path, str]): Location of the split file. This can be:
        
            - A name of the split file, which is stored as ``ppiref.definitions.PPIREF_SPLITS_DIR / f'{location}.json'``.
            - A Path object representing the path to the JSON split file.
            - An absolute path string starting with ``'/'`` representing the path to the JSON split file.

    Returns:
        Path: Path object representing the path to the JSON split file.
    """
    if not isinstance(location, str) and not isinstance(location, Path):
        raise ValueError('Wrong `location` type.')
    elif isinstance(location, str):
        if location.startswith('/'):
            if not location.endswith('.json'):
                raise ValueError(f'Wrong `location` extension in {location}.')
            location = Path(location)
        else:
            location = (PPIREF_SPLITS_DIR / location).with_suffix('.json')
    else:  # isinstance(location, Path) 
        if location.suffix != '.json':
            raise ValueError(f'Wrong `location` extension {location.suffix}.')
    return location
