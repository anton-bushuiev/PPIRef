"""Miscellaneous utility functions."""
import string
import random
import os
import requests
import warnings
import zipfile
from pathlib import Path
from typing import Union

from tqdm import tqdm

from ppiref.definitions import PPIREF_DATA_DIR


def random_id(length: int = 20) -> str:
    """Generate a random ID of a given length."""
    return ''.join(random.choices(string.digits, k=length))


def list_to_chunks(lst: list, n: int):
    """Cut list ``lst`` into ``n`` chunks."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_partition(lst: list, beg: float, end: float):
    """Get [``beg``, ``end``) partition of the list ``lst``."""
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
    project_url: str = 'https://zenodo.org/records/12821413/files/',
    destination_folder: Union[Path, str] = None
) -> None:
    """
    Download a file from Zenodo and extract its contents.

    Args:
        file (str): Name of the file to download and unpack. For example, ``'ppi_6A.zip'`` to
            download all the 6A-interfaces (based on the 6 Angstrom cutoff).
        project_url (str, optional): URL of the Zenodo project.
        destination_folder (Union[Path, str], optional): Path to the destination folder. If None, 
            the folder will be created in the ``ppiref.definitions.PPIREF_DATA_DIR`` directory. 
            Defaults to None.
    """
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

    # Download the file with progress bar
    response = requests.get(url, stream=True)
    total_size_in_bytes= int(response.headers.get('content-length', 0))
    block_size = 1024  # 1 Kibibyte
    progress_bar = tqdm(total=total_size_in_bytes, desc=f'Downloading to {destination_folder}', unit='iB', unit_scale=True)
    file_path = os.path.join(destination_folder, f'{stem}.zip')
    with open(file_path, 'wb') as file:
        for data in response.iter_content(block_size):
            progress_bar.update(len(data))
            file.write(data)
    progress_bar.close()

    # Check if the download was successful
    if total_size_in_bytes != 0 and progress_bar.n != total_size_in_bytes:
        raise RuntimeError("Download failed. Please try to download the file manually.")

    # Extract the contents of the zip file with progress bar
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        # List of archive contents
        list_of_files = zip_ref.infolist()
        total_files = len(list_of_files)
        # Progress bar for extraction
        with tqdm(total=total_files, desc="Extracting", unit='files') as pbar:
            for file in list_of_files:
                zip_ref.extract(member=file, path=destination_folder)
                pbar.update(1)

    # Delete the zip file after extraction
    os.remove(file_path)
