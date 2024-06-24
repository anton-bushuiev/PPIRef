"""Module for retrieving similar protein-protein interactions from large datasets."""
import subprocess
import tempfile
from pathlib import Path
from typing import Union, Optional

import pandas as pd

from ppiref.split import read_fold
from ppiref.definitions import PPIREF_DATA_DIR


class MMSeqs2PPIRetriever:
    def __init__(
        self,
        db: Optional[Union[str, Path]] = None,
        ppi_split: str = 'ppiref_6A_filtered',
        ppi_fold: str = 'whole',
        verbose: bool = False,
    ):
        """Retriever of similar protein-protein interactions using MMseqs2 sequence similarity
        search. 
        
        This class retrieves PPIs from a larger database containing similar sequences to any of the
        sequences involved in the query PPI.

        Please follow the official `MMseqs2 installation guide <https://github.com/soedinglab/mmseqs2?tab=readme-ov-file#installation>`_
        to install MMseqs2, and cite the original paper if you find the wrapper useful:

        .. code-block:: bibtex

            @article{steinegger2017mmseqs2,
                title={MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets},
                author={Steinegger, Martin and S{\"o}ding, Johannes},
                journal={Nature biotechnology},
                volume={35},
                number={11},
                pages={1026--1028},
                year={2017},
                publisher={Nature Publishing Group US New York}
            }

        Args:
            db (Union[str, Path], optional): MMseqs2 data base. Please see the MMseqs2 documentation
                on how to create a data base. Defaults to ``PPIRef/ppiref/ppi_6A_stats/mmseqs_db/db``,
                which indexes all protein sequences present in the PPIRef50K (6A version) dataset,
                compirising all proper protein-protein interactions from PDB.
            ppi_split (str, optional): Split of PPI ids to consider. In combination with the 
                ``ppi_fold`` argument, specifies the subset of PPIs to consider from the data base.
                Defaults to ``'ppiref_6A_filtered'``, which corresponds to complete PPIRef50K.
            ppi_fold (str, optional): Subset (fold) of PPI ids to consider from the specified split. 
                Defaults to ``'whole'`` to consider all PPIs from the split.
            verbose (bool, optional): If set to True prints MMseqs2 log. Defaults to False.
        """
        if db is None:
            db = PPIREF_DATA_DIR / 'ppiref/ppi_6A_stats/mmseqs_db/db'
        self.db = Path(db)
        assert db.is_file(), f'MMSeqs2 databsase file {db} does not exist.'
        self.verbose = verbose

        # Store all PPI ids in tuples (ppi_id, pdb_id, partners)
        ppis = read_fold(ppi_split, ppi_fold, full_paths=False)
        self.ppis = [(ppi, ppi.split('_', 1)[0], ppi.split('_')[1:]) for ppi in ppis]

    def query(self, seq: Union[str, Path]) -> tuple[list[float], list[str], list[str]]:
        """Query the MMseqs2 data base for protein-protein interactions with similar sequences.

        Args:
            seq (Union[str, Path]): Path to a query sequence in the fasta format.

        Returns:
            tuple[list[float], list[str], list[str]]: Tuple with three lists storing 3-tuples of
            PPI sequence similarity matches. The first list contains MMseqs2 sequence similarity 
            scores, the second list contains ids of the corresponding matched PPI entries, and 
            the third list contains chain names of the matched sequences from the corresponding 
            PPI entries.
        """
        # Prepare query sequence
        if isinstance(seq, str):
            if seq.endswith('.fasta'):
                seq = Path(seq)
            else:                
                raise NotImplementedError('Only fasta files are currently supported.')
        seq = Path(seq)
        assert seq.is_file(), f'File {seq} does not exist.'

        # Run MMSeqs2
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Prepare mmseqs2 arguments
            tmp_dir = Path(tmp_dir)
            result_file = tmp_dir / 'result.m8'

            # Run
            command = f'mmseqs easy-search {seq} {self.db} {result_file} {tmp_dir}'
            subprocess.run(command, check=False, capture_output=not self.verbose, shell=True)
            result = pd.read_csv(result_file, sep='\t', header=None)
            result = result.sort_values(by=2, ascending=False)

        # Map sequence similarity matches to all PPI similarity matches
        ps = result.loc[:, 1].to_list()
        seq_ids = result.loc[:, 2].to_list()
        ppis, ppi_seq_ids, partners = [], [], []
        for p, seq_id in zip(ps, seq_ids):
            pdb, partner = p.split(':')
            pdb = pdb.lower()
            for ppi_src, pdb_src, partners_src in self.ppis:
                if pdb == pdb_src and partner in partners_src:
                    ppis.append(ppi_src)
                    ppi_seq_ids.append(seq_id)
                    partners.append(partner)

        return ppi_seq_ids, ppis, partners
