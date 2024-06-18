import os
import io
import shutil
import warnings
import itertools
import random
import traceback
import multiprocessing
import subprocess
import concurrent.futures
import tempfile
from abc import ABC
from typing import Iterable, Callable, Literal, Optional, Union
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from Bio import Align
from Bio.Align import substitution_matrices
from tqdm import tqdm
from graphein.protein.config import ProteinGraphConfig
from graphein.protein.graphs import construct_graph
from graphein.protein.edges.distance import add_k_nn_edges
from graphein.protein.features.nodes.amino_acid import (
    amino_acid_one_hot, meiler_embedding)
from graphein.protein.features.sequence.embeddings import esm_residue_embedding
from graphein.protein.features.sequence.utils import (
    aggregate_feature_over_chains, aggregate_feature_over_residues)

from mutils.pdb import get_sequences

from ppiref.split import read_fold
from ppiref.utils.pdb import download_pdb
from ppiref.utils.ppipath import path_to_pdb_id, ppi_id_to_nested_suffix, path_to_partners
from ppiref.definitions import IALIGN_PATH, USALIGN_PATH, PPIREF_DATA_DIR


class PPIComparator(ABC):
    def __init__(self,
        max_workers: int = os.cpu_count() - 2,
        parallel_kind: Literal['threads', 'processes'] = 'processes',
        verbose=False
    ):
        """Abstract class for comparing protein-protein interactions (PPIs).

        Args:
            max_workers (int, optional): Number of workers to use for parallel operations (such as
                 comparing large sets of PPIs pairwise). Defaults to os.cpu_count() - 2.
            parallel_kind (Literal[&#39;threads&#39;, &#39;processes&#39;], optional): Use
                 multi-treading or multi-processing for parallel operations. Defaults to 'processes'.
            verbose (bool, optional): If set to True, prints detailed log to the standard output.
                 May be useful for debugging. Defaults to False.
        """
        self.max_workers = max_workers
        self.parallel_kind = parallel_kind
        self.verbose = verbose

    def compare(self, ppi0: Path, ppi1: Path) -> dict:
        """Abstract method for comparing two PPIs. Should be implemented in a subclass.

        Args:
            ppi0 (Path): Path to the first .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.
            ppi1 (Path): Path to the second .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.

        Returns:
            dict: Dictionary with comparison results.
        """
        raise NotImplementedError()

    def compare_all_against_all(
        self,
        ppis0: Iterable[Path] = None,
        ppis1: Iterable[Path] = None,
        ppi_pairs: Iterable[Iterable[Path]] = None
    ) -> pd.DataFrame:
        """Comparing all PPIs from one set against all PPIs from another set. This method from the
         abstract class is used as a default implementation where all PPI pairs are compared 
         in a data parallel way. A subclass may override this method to provide a more efficient
         implementation.

        Args:
            ppis0 (Iterable[Path], optional): First set of PPI paths. Defaults to None.
            ppis1 (Iterable[Path], optional): Second set of PPI paths. Defaults to None.
            ppi_pairs (Iterable[Iterable[Path]], optional): Pre-defined pairs to compare instead of
                 complete pair-wise comparison of two sets. Defaults to None.

        Returns:
            pd.DataFrame: Data frame with comparison results. The data frame has two columns
                 corresponding to pairs of PPI ids, and additional columns with comparison metrics.
        """
        # Parse input
        if ppis0 is not None and ppis1 is not None:
            ppi_pairs = list(itertools.product(ppis0, ppis1))
        elif ppi_pairs is None:
            raise ValueError('Input pairs are not specified.')

        # Compare and store to dataframe
        df = self._execute_task_parallel(
            self.compare, ppi_pairs, desc=f'Comparing PPIs with {self.__class__.__name__}'
        )
        df = pd.DataFrame(df)
        return df

    def _execute_task_parallel(
        self,
        func: Callable,  
        inputs: Iterable,
        kind: str = None,
        desc: str = '',
        chunksize: int = 1000
    ) -> Iterable:
        """Universal method for executing a function in parallel on a set of inputs. The method
         uses either multi-threading or multi-processing depending on the `parallel_kind` attribute
         of the class.

        Args:
            func (Callable): Function to apply to inputs
            inputs (Iterable): All inputs to apply the function to
            kind (str, optional): Use multi-threading or multi-processing. Defaults to None to use
                 the kind specified in the comparator constructor.
            desc (str, optional): Progress bar description. Defaults to ''.
            chunksize (int, optional): Number of inputs to be processed at a time by a single
                 worker. Low chunksize may result into a too large number (e.g. >100K) of jobs for
                 workers, which may lead to a significant management overhead. Defaults to 1000.

        Returns:
            Iterable: List of results of applying the function to inputs.
        """
        # Init executor
        if kind is None:
            kind = self.parallel_kind
        if kind == 'threads':
            executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers)
        elif kind == 'processes':
            executor = concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers)
        else:
            raise ValueError("Invalid 'kind'. Use 'threads' or 'processes'.")

        # Run tasks
        total_tasks = len(inputs)
        futures = []
        with tqdm(desc=f'{desc} ({executor._max_workers} {kind})', total=total_tasks) as pbar:
            # Submit tasks
            for inp in inputs:
                future = executor.submit(self._unpacked_call, func, inp)
                futures.append(future)
            
            results = []
            # Update progress bar as tasks complete
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
                pbar.update(1)
        
        return results
    
    def _unpacked_call(self, func, args):
        """Helper function to unpack arguments for a function call. Used for parallel execution."""
        return func(*args)
    

class USalign(PPIComparator):
    def __init__(
        self,
        path: Path = USALIGN_PATH,
        args: str = '',
        **kwargs
    ):
        """Wrapper for the US-align protein-protein interaction comparator
         (https://zhanggroup.org/US-align/). Please install the US-align software and provide the
         path to the executable.

        Args:
            path (Path, optional): Path to the `USalign` executable. Defaults to USALIGN_PATH, which
                 should be used if you follow the instructions from the PPIRef/external/README.md to
                 install US-align.
            args (str, optional): Optional command line arguments to be passed to US-align.
                 Defaults to ''.
        """
        super().__init__(**kwargs)
        self.path = path
        self.args = args

    def compare(self, ppi0: Path, ppi1: Path) -> dict:
        """Compare two protein-protein interactions.

        Args:
            ppi0 (Path): Path to the first .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.
            ppi1 (Path): Path to the second .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.

        Returns:
            dict: Dictionary with two PPI ids being compared and all comparison metrics produced
                 by US-align.
        """
        ppi_id0, ppi_id1 = ppi0.stem, ppi1.stem

        command = f'{self.path} {ppi0} {ppi1} -outfmt 2 -mm 1 -ter 1 {self.args}'
        out = subprocess.check_output(command, capture_output=False, shell=True)
        out = out.decode()
        metrics = pd.read_csv(io.StringIO(out), sep='\t')
        metrics = metrics.to_dict('records')[0]
        del metrics['#PDBchain1']
        del metrics['PDBchain2']

        return {'PPI0': ppi_id0, 'PPI1': ppi_id1} | metrics


class IAlign(PPIComparator):
    def __init__(
        self,
        path: Path = IALIGN_PATH,
        args: str = '-a 2 -minp 1 -mini 1',
        tmp_out_pref: str = 'iAlign_tmp_out',
        multiple_resolution: tuple[str, Literal['max', 'min']] = ('P-value', 'min'),
        **kwargs
    ):
        """Wrapper for the iAlign protein-protein interaction comparator
         (https://doi.org/10.1093/bioinformatics/btq404). iAlign is the adaptation of TM-align to
         protein-protein interfaces. Please install the iAlign source code in Perl and provide the
         path to the executable.

        Args:
            path (Path, optional): Path to the `ialign.pl` executable. Defaults to USALIGN_PATH, which
                 should be used if you follow the instructions from the PPIRef README to install
                 iAlign.
            args (str, optional): Command line arguments to use for iAlign. Defaults to
                 '-a 2 -minp 1 -mini 1' to operate in verbose mode and not to skip small interfaces.
            tmp_out_pref (str, optional): Prefix to use for temporariry directories storing iAlign
                 outputs. Defaults to 'iAlign_tmp_out'.
            multiple_resolution (tuple[str, Literal['max', 'min']], optional): If one or both input
                 files contain multiple interfaces, this parameter specifies how to resolve multiple
                 comparison results. The first element of the tuple specifies the metric to use for
                 resolution, and the second element specifies whether to use the maximum or minimum.

        Notes:
            TODO: Remove `tmp_out_pref` and use the `tempdir` package instead.
        """
        super().__init__(**kwargs)
        self.path = path
        self.args = args
        self.tmp_out_pref = tmp_out_pref
        self.multiple_resolution = multiple_resolution

    def compare(self, ppi0: Path, ppi1: Path) -> dict:
        """Compare two protein-protein interactions.

        Args:
            ppi0 (Path): Path to the first .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.
            ppi1 (Path): Path to the second .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.

        Returns:
            dict: Dictionary with two PPI ids being compared and all comparison metrics produced
                 by iAlign.
        """
        # Create temporary directory with unique name to store outputs
        ppi_id0, ppi_id1 = ppi0.stem, ppi1.stem
        tmp_out_dir = Path(f'{self.tmp_out_pref}_{ppi_id0}_{ppi_id1}')
        tmp_out_dir.mkdir(parents=True, exist_ok=True)
        tmp_out_file = tmp_out_dir / Path(f'{ppi_id0}_{ppi_id1}.out')

        # Create a temporary copy of one of the files if necessary
        # (iAlign does not handle comparison of files with same pdb stem)
        same_pdb = False
        ppi1_old = ppi1
        if path_to_pdb_id(ppi0) == path_to_pdb_id(ppi1):
            same_pdb = True
            pdb_id, pair_id = ppi1.stem.split('_', 1)
            tmp_stem = f'{pdb_id}-{random.randint(1000, 9999)}.{pair_id}'
            tmp_name = ppi1.with_stem(tmp_stem).name
            ppi1 = tmp_out_dir / tmp_name
            shutil.copy(ppi1_old, ppi1)

        # Run iAlign
        pdb_id0_parsed, pdb_id1_parsed = ppi_id0.split('_'), ppi_id1.split('_')
        if len(pdb_id0_parsed) > 1 and len(pdb_id1_parsed) > 1:
            chains = f"-c1 {''.join(pdb_id0_parsed[1:])} -c2 {''.join(pdb_id1_parsed[1:])}"
        else:
            chains = ''
        command = (
            f'perl {self.path} {self.args} -w {tmp_out_dir}'
            f' {chains}'
            f' {ppi0} {ppi1} > {tmp_out_file}'
        )
        if self.verbose:
            print(command)
        subprocess.run(command, check=False, capture_output=not self.verbose, shell=True)

        # Parse output
        with open(tmp_out_file, 'r') as file:
            out_lines = list(map(lambda x: x.strip(), file.readlines()))

        # Parse all metric segments
        # Example:
        # >>>3W2D_A_H_LHL vs 4Y61_A_BAB
        #
        # Structure 1: 3W2D_A_H_LHL_int.pdb Chains H L,  83 AAs,  95 Contacts
        # Structure 2: 4Y61_A_BAB_int.pdb   Chains A B,  35 AAs,  36 Contacts
        #
        # IS-score = 0.14796, P-value = 0.5985E+000, Z-score =  0.092
        # Number of aligned residues  =  14
        # Number of aligned contacts  =   4
        # RMSD =  3.43, Seq identity  = 0.071
        start_metric_lines = [i for i, line in enumerate(out_lines) if line.startswith('>>>')]
        pairwise_metrics = []
        for start_metric_line in start_metric_lines:
            metric_lines = list(range(start_metric_line + 5, start_metric_line + 5 + 4))
            metric_lines = list(zip(metric_lines, [float, int, int, float]))

            metrics = {}
            if metric_lines[-1][0] >= len(out_lines):
                if self.verbose:
                    print(
                        f'iAlign error comparing ({ppi0}, {ppi1_old})\n'
                        f'iAlign output:\n' + '\n'.join(out_lines),
                        flush=True
                    )
                    break
            else:
                for line_id, dtypes in metric_lines:
                    line = out_lines[line_id]
                    metrics |= self._parse_metric_line(line, dtypes)
                if self.verbose:
                    print(f'Finished comparing ({ppi0}, {ppi1_old})', flush=True)
            pairwise_metrics.append(metrics)
        
        # Select the comparison with the best hit
        if pairwise_metrics:
            resolution_func = max if self.multiple_resolution[1] == 'max' else min
            metrics = resolution_func(pairwise_metrics, key=lambda x: x[self.multiple_resolution[0]])
        else:
            metrics = {}

        # Clean
        if same_pdb:
            ppi1.unlink()
        tmp_out_file.unlink()
        shutil.rmtree(tmp_out_dir)

        # Return result dict
        return {'PPI0': ppi_id0, 'PPI1': ppi_id1} | metrics

    @staticmethod
    def _parse_metric_line(line: str, dtypes=float) -> dict:
        """Parse a line with a metric from iAlign output."""
        metrics = map(lambda x: x.split('='), line.split(','))
        metrics = {metric.strip(): dtypes(val) for metric, val in metrics}
        return metrics


class SequenceIdentityComparator(PPIComparator):
    def __init__(
        self,
        pdb_dir: Union[Path, str],
        nested_pdb_dir: bool = False,
        aligner: Align.PairwiseAligner = None,
        **kwargs
    ):
        """Protein-protein interaction comparator based on sequence identity. The comparator uses 
         the BioPython library to align protein sequences and calculate the pairwise sequence
         identity. The similarity is calculated as the maximum pairwise sequence identity between
         chains from different PPIs.

        Args:
            pdb_dir (Union[Path, str]): Directory with complete .pdb files that were used to extract
                 PPI interactions from. The directory is used to get complete protein sequences for
                 comparison.
            nested_pdb_dir (bool): True if files are in the format `pdb_dir/bc/abcd.pdb`.
                 False if files are in the format `pdb_dir/abcd.pdb`. Defaults to False.
            aligner (Align.PairwiseAligner, optional): BioPython Aligner to use. Defaults to None
                 to use the one employed in the PoseBusters package (https://arxiv.org/abs/2308.05777).
        """
        super().__init__(**kwargs)
        self.pdb_dir = Path(pdb_dir)
        self.nested_pdb_dir = nested_pdb_dir

        if aligner is None:
            aligner = Align.PairwiseAligner()
            aligner.open_gap_score = -11
            aligner.extend_gap_score = -1
            aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        self.aligner = aligner

    def compare(self, ppi0: Path, ppi1: Path) -> dict:
        """Compare two protein-protein interactions.

        Args:
            ppi0 (Path): Path to the first .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.
            ppi1 (Path): Path to the second .pdb file containing a PPI. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.

        Returns:
            dict: Dictionary with two PPI ids being compared and maximum pairwise sequence identity
                 between chains from different PPIs.
        """
        # Parse input file names
        ppi_id0, ppi_id1 = ppi0.stem, ppi1.stem
        pdb_id0, pdb_id1 = path_to_pdb_id(ppi0), path_to_pdb_id(ppi1)
        partners0 = path_to_partners(ppi0)
        partners1 = path_to_partners(ppi1)

        # Get paths to complete .pdb files (with full chains)
        if self.nested_pdb_dir:
            pdb0 = self.pdb_dir / ppi_id_to_nested_suffix(pdb_id0)
            pdb1 = self.pdb_dir / ppi_id_to_nested_suffix(pdb_id1)
        else:
            pdb0 = self.pdb_dir / (path_to_pdb_id(ppi0) + '.pdb')
            pdb1 = self.pdb_dir / (path_to_pdb_id(ppi1) + '.pdb')

        # Get full sequences for each PPI
        seqs0 = get_sequences(pdb0)
        seqs1 = get_sequences(pdb1)

        # Calculate pairwise sequence identity for all pairs of chains from different PPIs
        pairwise = []
        for p0 in partners0:
            for p1 in partners1:
                seq_id = self._sequence_identity(seqs0[p0], seqs1[p1])
                pairwise.append(seq_id)

        metrics = {'Maximum pairwise sequence identity': max(pairwise)}
        return {'PPI0': ppi_id0, 'PPI1': ppi_id1} | metrics
    
    def _sequence_identity(self, seq0: str, seq1: str) -> float:
        """Calculate sequence identity between two protein sequences with a BioPython Aligner."""
        alignments = self.aligner.align(seq0, seq1)
        alignment = next(alignments)
        
        parts = alignment.format("fasta").split('\n')
        seq0, seq1 = parts[1], parts[3]
        matches = sum(1 for a, b in zip(seq0, seq1) if a == b)
        percent_match = matches / min(len(seq0), len(seq1))
        return percent_match


IDIST_EMBEDDING_KIND = Literal[
    'amino_acid_one_hot', 'esm_embedding', 'meiler_embedding'
]


class IDist(PPIComparator):
    MAX_INTERFACE_SIZE = 1_000_000

    def __init__(
        self,
        kind: IDIST_EMBEDDING_KIND = 'amino_acid_one_hot',
        near_duplicate_threshold: float = 0.04,
        pdb_dir: Optional[Union[str, Path]] = None,
        *args,
        **kwargs
    ):
        """Implementation of iDist protein-protein interaction comparator
        (https://arxiv.org/abs/2310.18515).
        
        The comparator uses a simple non-parametrized one-step
         message passing to embed protein-protein interfaces and then compares them using Euclidean
         distance. iDist approximates 3D alignment-based methods, iAlign and US-align, on detecting
         near-duplicate protein-protein interfaces. iDist is more than 100 times faster and finds
         same near duplicates with 99% precision and 97% recall.

        Args:
            kind (IDIST_EMBEDDING_KIND, optional): Kind of node embeddings to use for message
                 passing. Defaults to 'amino_acid_one_hot' which leads to the best alignment 
                 approximation performance.
            near_duplicate_threshold (float, optional): Threshold on Euclidean distance to detect
                 near-duplicate interfaces. It is recommended to use the threshold of 0.04 for the
                 interfaces extracted with the 6A cutoff radius between heavy atoms, and the 
                 threshold of 0.03 with the 10A cutoff. Please see the paper for details. Defaults
                 to 0.04.
            pdb_dir (Optional[Path], optional): Directory storing complete .pdb files that were used
                 to exctract interfaces from. Should be not None if `kind == 'esm_embedding'`, as the
                 ESM protein language model is used with full protein sequences. Defaults to None.
        """
        super().__init__(*args, **kwargs)

        # Prepare node features construction
        self.kind = kind
        if kind == 'amino_acid_one_hot':
            dimer_node_metadata_functions = []
            ppi_node_metadata_functions = [amino_acid_one_hot]
        elif kind == 'meiler_embedding':
            dimer_node_metadata_functions = []
            ppi_node_metadata_functions = [meiler_embedding]
        elif kind == 'esm_embedding':
            dimer_node_metadata_functions = []
            ppi_node_metadata_functions = []
        else:
            raise ValueError('Unknown `kind` value.')

        # Init Graphein graph construction config for interfaces (or PPIs)
        self.graphein_ppi_config = ProteinGraphConfig(
            edge_construction_functions=[
                partial(add_k_nn_edges, k=self.MAX_INTERFACE_SIZE,
                        long_interaction_threshold=0,
                        exclude_edges=['inter'], kind_name='intra'),
                partial(add_k_nn_edges, k=self.MAX_INTERFACE_SIZE,
                        long_interaction_threshold=0,
                        exclude_edges=['intra'], kind_name='inter')
            ],
            node_metadata_functions=ppi_node_metadata_functions,
            insertions=True
        )

        # Init Graphein graph construction config for complete complexes (or dimers)
        self.graphein_dimer_config = ProteinGraphConfig(
            edge_construction_functions=[],
            node_metadata_functions=dimer_node_metadata_functions,
            insertions=True,
        )

        # Init near-duplicate threshold
        self.near_duplicate_threshold = near_duplicate_threshold

        # Init directory for complete (or dimer) PDB files
        self.pdb_dir = pdb_dir
        if self.pdb_dir is not None:
            self.pdb_dir = Path(self.pdb_dir)
            if not self.pdb_dir.exists():
                self.pdb_dir.mkdir(parents=True, exist_ok=True)

        # Init cache for embeddings
        self.embeddings = dict()

        # Init index for near-duplicate detection
        self.neigh = None

    def compare(
        self,
        path0: Union[Path, str],
        path1: Union[Path, str]
    ) -> dict:
        """Compare two protein-protein interfaces.

        Args:
            ppi0 (Union[Path, str]): Path to the first .pdb file containing a PPI. It is recommended
                 to use the files produced by the `ppiref.extraction.PPIExtractor` class.
            ppi1 (Union[Path, str]): Path to the second .pdb file containing a PPI. It is 
                 recommended to use the files produced by the `ppiref.extraction.PPIExtractor`
                 class.

        Returns:
            dict: Dictionary with two PPI ids being compared and the resulting iDist distance
                 (Euclidean distance in the embedding space).
        """
        path0, path1 = Path(path0), Path(path1)
        pdb0, pdb1 = path0.stem, path1.stem

        # Encode and compare
        emb0 = self.embed(path0)
        emb1 = self.embed(path1)
        metrics = {
            'iDist': np.linalg.norm(emb0 - emb1),
            # 'L1': np.linalg.norm(emb0 - emb1, ord = 1),
            # 'Cosine Similarity':
            #     np.dot(emb0, emb1) / (np.linalg.norm(emb0) * np.linalg.norm(emb1))
        }

        # Return result dict
        return {'PPI0': pdb0, 'PPI1': pdb1} | metrics

    # TODO Accelerate with sklearn pairwise
    def compare_all_against_all(
        self,
        ppis0: Iterable[Path] = None,
        ppis1: Iterable[Path] = None,
        ppi_pairs: Iterable[Iterable[Path]] = None,
        embed: bool = True
    ) -> pd.DataFrame:
        """Compare all PPIs from one set against all PPIs from another set efficiently using iDist.

        Args:
            ppis0 (Iterable[Path], optional): First set of PPI paths. Defaults to None.
            ppis1 (Iterable[Path], optional): Second set of PPI paths. Defaults to None.
            ppi_pairs (Iterable[Iterable[Path]], optional): Pre-defined pairs to compare instead of
                 complete pair-wise comparison of two sets. Defaults to None.
            embed (bool, optional): If set to True, embeds all PPIs before comparison not to repeat
                 same embeddign twice. Defaults to True.

        Returns:
            pd.DataFrame: Data frame with comparison results. The data frame has two columns
                 corresponding to pairs of PPI ids, and an additional column with iDist distances.
        """
        # Parse input
        if ppis0 is not None and ppis1 is not None:
            ppi_pairs = list(itertools.product(ppis0, ppis1))
        elif ppi_pairs is None:
            raise ValueError('Input pairs are not specified.')

        # Embed PPIs
        if embed:
            if ppis0 is not None and ppis1 is not None:
                ppis = set(ppis0) | set(ppis1)
            else:  # ppi_pairs is not None
                ppis = set(itertools.chain(*ppi_pairs))
            self.embed_parallel(ppis)

        # Compare PPIs pairwise
        df = [self.compare(*x) for x in ppi_pairs]
        df = pd.DataFrame(df)
        return df

    def embed(self, ppi: Path, store: bool = True) -> np.array:
        """Embed a protein-protein interface.

        Args:
            ppi (Path): Path to the .pdb file containing a PPI. It is recommended to use the files
                 produced by the `ppiref.extraction.PPIExtractor` class.
            store (bool, optional): Set to True to store the embedding in cache. Defaults to True.

        Returns:
            np.array: PPI interface embedding.
        """
        ppi_id = ppi.stem
        ppi = str(ppi)
        if self.kind == 'esm_embedding':
            dimer = self._ppi_to_pdb(ppi)
            dimer = str(dimer)

        # Get embedding from cache
        if ppi_id in self.embeddings:
            return self.embeddings[ppi_id]

        # Construct interface graph
        g_ppi = construct_graph(
            config=self.graphein_ppi_config, path=ppi, verbose=False
        )

        # Construct complete dimer graph and transfer embeddings to interface subgraph if ESM
        # embeddings are used
        if self.kind == 'esm_embedding':
            # Construct protein graph for complete dimer structure
            g_dimer = construct_graph(
                config=self.graphein_dimer_config, pdb_path=dimer,
                chain_selection=''.join(g_ppi.graph['chain_ids']), verbose=False
            )
            # Add ESM residue embeddings to complete-structure graph
            g_dimer = esm_residue_embedding(g_dimer)

            # Transfer node embeddings from complete structure to PPI subgraph
            for n in g_ppi.nodes():
                g_ppi.nodes[n]['esm_embedding'] = \
                    g_dimer.nodes[n]['esm_embedding']
                
        # Rename node features if Meiler embeddings are used
        elif self.kind == 'meiler_embedding':
            for v in g_ppi.nodes():
                g_ppi.nodes[v]['meiler_embedding'] = \
                    g_ppi.nodes[v]['meiler'].values

        # Aggregate neighborhoods (non-parametrized one-step message passing)
        # Note: Can be significantly accelerated via graph-level matmuls
        for v in g_ppi.nodes():
            msg_inter = []
            msg_intra = []
            for n, e in g_ppi[v].items():
                signal = \
                    np.exp(-(e['distance']/4)**2) * g_ppi.nodes[n][self.kind]
                if 'inter' in e['kind']:
                    msg_inter.append(-signal)
                elif 'intra' in e['kind']:
                    msg_intra.append(signal)
            msg_inter = np.mean(msg_inter, axis=0)
            msg_intra = np.mean(msg_intra, axis=0)
            msg = np.mean([msg_inter, msg_intra], axis=0)
            g_ppi.nodes[v]['embedding'] = np.mean([
                g_ppi.nodes[v][self.kind],
                msg
            ], axis=0)

        # Aggregate residues in chains and then chain embeddings
        aggregate_feature_over_residues(g_ppi, 'embedding', 'mean')
        aggregate_feature_over_chains(g_ppi, 'embedding_mean', 'mean')

        # Save to cache and return
        embedding = g_ppi.graph['embedding_mean_mean']
        if store:
            self.embeddings[ppi_id] = embedding
        return embedding

    def embed_without_exception(self, ppi: Path) -> np.array:
        """Embed a PPI and catch exceptions to avoid breaking the parallel execution.

        Args:
            ppi (Path): Path to the .pdb file containing a PPI. It is recommended to use the files
                 produced by the `ppiref.extraction.PPIExtractor` class.

        Returns:
            np.array: PPI interface embedding.
        """
        try:
            embedding = self.embed(ppi)
        except Exception as exc:
            print(f'{ppi} led to an exception {exc}:')
            print(traceback.format_exc(), end='\n\n')
            embedding = np.full(1024, np.nan)
        if self.verbose:
            print(f'Embedded {ppi}')
        return embedding

    def embed_parallel(self, ppis: Iterable[Path]) -> None:
        """Embed a set of PPIs in parallel and store in cache.

        Args:
            ppis (Iterable[Path]): Paths to the .pdb files containing PPIs. It is recommended to use
                 the files produced by the `ppiref.extraction.PPIExtractor` class.
        """
        # Adapt dict for multi-processing
        if self.parallel_kind == 'processes':
            self.embeddings = multiprocessing.Manager().dict()

        # Embed in parallel
        ppis = list(map(lambda x: (x,), ppis))
        self._execute_task_parallel(
            self.embed_without_exception, ppis, desc='Embedding PPIs'
        )

        # Return dict back to ordinary
        self.embeddings = dict(self.embeddings)

    def deduplicate_embeddings(self) -> None:
        """Deduplicate embeddings based on the threshold Euclidean distance between them. 
        
        The method iteratively removes embeddings that are closer than the threshold to any other
         embedding while making only one-sided comparisons (i.e. a<->b but not b<->a). Since the 
         number of embeddings may be large, the method processes the pairwise distances matrix
         in chunks of consecuite rows. 
        """
        df_emb = self.get_embeddings()
        pad_val = -1

        # Process adjacency chunk and return duplicated ids
        def reduce_func(chunk, start):
            chunk = chunk < self.near_duplicate_threshold
            chunk &= ~np.tri(*chunk.shape, k=start).astype(bool)
            idx = chunk.sum(axis=1).nonzero()[0]
            idx += start
            idx = np.pad(idx, (0, len(chunk) - len(idx)), constant_values=pad_val)
            return idx
    
        # Iterate over chunks of adjacency matrix
        def get_chunks():
            chunks = sklearn.metrics.pairwise_distances_chunked(
                df_emb,
                n_jobs=self.max_workers,
                working_memory=sklearn.get_config()['working_memory'],
                reduce_func=reduce_func
            )
            return chunks
        
        # Get chunk size
        chunk_size = len(next(get_chunks()))
        n_chunks = int(np.ceil(len(df_emb) / chunk_size))
        
        # Run deduplication over chunks
        idx_to_remove = []
        for chunk in tqdm(get_chunks(), total=n_chunks, desc='Processing adjacency chunks'):
            chunk = chunk[chunk != pad_val]
            idx_to_remove.extend(chunk)
        names_to_remove = df_emb.index[idx_to_remove]

        # Convert to original dict format
        self.embeddings = {
            name: z for name, z in self.embeddings.items() if name not in names_to_remove
        }

    def build_index(self) -> None:
        """Build an index for fast near-duplicate detection based on Euclidean distance between 
        embeddings.
        """
        self.neigh = sklearn.neighbors.NearestNeighbors(radius=self.near_duplicate_threshold)
        self.neigh.fit(self.get_embeddings())

    def query(self, q: np.array) -> np.array:
        """Query the index for near-duplicate embeddings.

        Args:
            q (np.array): Input query embedding(s). May be a single embedding or a stack of multiple
                 embeddings.

        Returns:
            np.array: If a single query embedding is provided, returns an array of near-duplicate
                 PPI ids. If multiple query embeddings are provided, returns an array of arrays of
                 near-duplicate PPI ids. The ids are sorted by their distance to the query embedding(s),
                 with the first id being the closest near duplicate.
        """
        # Build index if not built yet
        if self.neigh is None:
            self.build_index()

        # Query index
        assert len(q.shape) in [1, 2]
        if len(q.shape) == 1:  # single query vector
            neigh_dist, neigh_ind = self.neigh.radius_neighbors(np.expand_dims(q, 0), sort_results=True)
            neigh_dist, neigh_ind = neigh_dist[0], neigh_ind[0]
            names = self.get_embeddings().index
            neigh_ind = names[neigh_ind].to_numpy()
        else:  # multiple stacked query vectors
            neigh_dist, neigh_ind = self.neigh.radius_neighbors(q, sort_results=True)
            names = self.get_embeddings().index
            neigh_ind = [names[row].to_numpy() for row in neigh_ind]
            
        return neigh_dist, neigh_ind

    def get_embeddings(self) -> pd.DataFrame:
        """Get the embeddings stored in the cache as a Pandas data frame.

        Returns:
            pd.DataFrame: Data frame with embeddings in rows indexed by PPI ids.
        """
        return pd.DataFrame(dict(self.embeddings)).T
    
    def write_embeddings(self, path: Path) -> None:
        """Write the embeddings stored in the cache to a .csv file."""
        self.get_embeddings().to_csv(path)

    def read_embeddings(self, df: Union[Path, pd.DataFrame], dropna: bool = False) -> None:
        """Read embeddings from a .csv file or a Pandas data frame and store them in the cache.

        Args:
            df (Union[Path, pd.DataFrame]): Pandas data frame with embeddings.
            dropna (bool, optional): Drop rows (embeddings) containing NaN values. NaN values may
                 appear if iDist fails to embed some PPIs. Defaults to False.
        """
        if isinstance(df, Path):
            df = pd.read_csv(df, index_col=0)
        if dropna:
            df = df.dropna()
        embeddings = df.T.to_dict(orient='series')
        embeddings = {k: np.array(v) for k, v in embeddings.items()}
        self.embeddings = embeddings

    def _ppi_to_pdb(self, ppi: Path) -> Path:
        """Convert a path to a PPI interface to a path to a complete PDB file based on the directory
         specified in the `pdb_dir` class attribute. If files are not found in the directory, they
         are downloaded from the PDB database.

        Args:
            ppi (Path): Path to the .pdb file containing a PPI interface. It is recommended to use 
                 the files produced by the `ppiref.extraction.PPIExtractor` class.

        Returns:
            Path: Path to the corresponding .pdb file containing the complete interaction structure.
        """
        if self.pdb_dir is None:
            raise ValueError('`pdb_dir` is not specified.')

        pdb_id = path_to_pdb_id(ppi)

        # pdb_dir has separated structure case (e.g. pdb_dir/cd/abcd.pdb)
        pdb = self.pdb_dir / f'{pdb_id[1:3]}' / f'{pdb_id}.pdb'
        if pdb.is_file():
            return pdb

        # pdb_dir has ordinary flattened structure case
        pdb = self.pdb_dir / f'{pdb_id}.pdb'
        if not pdb.is_file():
            download_pdb(pdb_id, path=pdb)

        return pdb


class MMSeqs2PPIRetriever:
    DEFAULT_DB = PPIREF_DATA_DIR / 'ppiref/ppi_6A_stats/mmseqs_db/db'

    def __init__(
        self,
        db: Union[str, Path] = DEFAULT_DB,
        ppi_split: str = 'ppiref_6A_filtered',
        ppi_fold: str = 'whole',
        verbose: bool = False,
    ):
        """Retriever of similar protein-protein interactions using MMseqs2 sequence similarity
             search (https://www.nature.com/articles/nbt.3988). Specifically, this class retrieves 
             PPIs containing similar sequences to any of the sequences involved in the query PPI.

        Args:
            db (Union[str, Path], optional): MMseqs2 data base. Please see the MMseqs2 documentation
                 on how to create a data base. Defaults to `MMSeqs2PPIRetriever.DEFAULT_DB`, which
                 indexes all protein sequences present in the PPIRef50K (6A version) dataset, 
                 compirising all proper protein-protein interactions from PDB.
            ppi_split (str, optional): Split of PPI ids to consider. In combination with the 
                 `ppi_fold` argument, specifies the subset of PPIs to consider from the data base.
                 Defaults to 'ppiref_6A_filtered', which corresponds to complete PPIRef50K.
            ppi_fold (str, optional): Subset (fold) of PPI ids to consider from the specified split. 
                 Defaults to 'whole' to consider all PPIs from the split.
            verbose (bool, optional): If set to True prints MMseqs2 log. Defaults to False.
        """
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
