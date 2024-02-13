import string
import random
import os
import requests
import warnings
import zipfile
from pathlib import Path
from typing import Union

from ppiref.definitions import PPIREF_DATA_DIR


def random_id(length: int = 20) -> str:
    return ''.join(random.choices(string.digits, k=length))


def list_to_chunks(lst: list, n: int):
    """Cut list `lst` into `n` chunks.
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_partition(lst: list, beg: float, end: float):
    """Get [`beg`, `end`) partition of the list `lst`.
    """
    if not 0. <= beg <= end <= 1.:
        raise ValueError(
            "Invalid range. beg must be less than or equal to end,"
            "and both must be between 0 and 1."
        )

    n = len(lst)
    start_index = int(n * beg)
    end_index = int(n * end)

    return lst[start_index:end_index]


def download_from_zenodo(
    file: str,
    project_url: str = 'https://zenodo.org/records/10651459/files/',
    destination_folder: Union[Path, str] = None
) -> None:
    # Create full file url
    url = project_url + file
    stem = Path(url).stem

    # Create fill path to the destination folder
    if destination_folder is None:
        destination_folder = PPIREF_DATA_DIR / 'ppiref' / stem

    # Check if the folder already exists
    if (destination_folder).is_dir():
        warnings.warn(f'{destination_folder} already exists. Skipping download.')
        return
    
    # Create the folder
    destination_folder.mkdir(parents=True, exist_ok=True)

    # Download the file
    response = requests.get(url)
    file_path = os.path.join(destination_folder, f'{stem}.zip')
    with open(file_path, 'wb') as file:
        file.write(response.content)

    # Extract the contents of the zip file
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(destination_folder)

    # Delete the zip file after extraction
    os.remove(file_path)
