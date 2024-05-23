<div align="center">

# PPIRef

</div>

<p align="center">
  <img src="assets/readme-dimer-close-up.png"/>
</p>

PPIRef is a Python package for working with 3D structures of protein-protein interactions (PPIs). It is based on the PPIRef dataset, comprising all PPIs from Protein Data Bank.

Please do not hesitate to contact us or create an issue/PR if you have any questions or suggestions ✌️.

# PPIRef dataset

<p align="center">
  <img align="right" width="350" src="https://github.com/anton-bushuiev/PPIRef/assets/67932762/f5d12ffb-8d1b-40b8-bab1-23d33a091f05"/>
</p>

PPIRef is a complete* and non-redundant dataset of protein-protein interactions (PPIs). It was constructed in three steps: (i) exhaustevily extract all putative protein dimers from PDB based on heavy atom contacts (we use two variants with 6A and 10A distance cutoffs), (ii) filter out not proper PPIs based on buried surface area and quality criteria, (iii) remove near-duplicate PPIs with iDist - a fast algorithm accurately approximating PPI alignment-based algorithms such as iAlign. See [our paper](https://arxiv.org/abs/2310.18515) for details.

\* (with respect to PDB in Jan, 2024)

## How to use

**Downloading**

The PPIRef dataset is available on [Zenodo](https://zenodo.org/records/10700674). We recommended downloading it and putting under `ppiref/data/ppiref`. This can be done automatically by running the following command to download all 6A interfaces in .pdb files (see below on how to install the ppiref pacakge):
```python
from ppiref.utils.misc import download_from_zenodo
download_from_zenodo('ppi_6A.zip')
```

**Reconstructing from scratch**

Alternatively, the dataset may be reconstructed from scratch with the following steps: (i) downloading and unpacking PDB (`scripts/download_pdb.sh` and `scripts/unpack_pdb.sh`), (ii) extracting all distance-based dimeric PPIs from all PDB files (`scripts/ppi_extractor.py`) and (iii) filtering and clustering PPIs using iDist (see the corresponding section below). The `ppi_6A.zip` archive contains the result of steps (i) and (ii), and the filtered subsets resulting from step (iii) are stored in `ppiref/data/splits`, as discussed below.

**Using subsets and splits**

After downloading or reconstructing the dataset, you will obtain all the dimeric interactions in PDB that have at least one contact between heavy atoms within the threshold distance (6 Angstrom for `ppi_6A.zip`). Then, you can use the predefined splits and other subsets, located in `ppiref/data/splits`. They can be read with the `ppiref.split.read_split` function (see the corresponding section below). The splits are stored in JSON files and contain the PPI IDs (e.g., `1bui_A_C` for the interaction between chains `A` and `C` in the `1bui` PDB entry) for each subset (e.g., `train`, `test`, `val`). Some examples are:
- `ppiref_6A_filtered.json` - all 6A interfaces that satisfy the filtering criteria for proper PPIs (a.k.a PPIRef300K, see Appendix A in the [paper](https://arxiv.org/abs/2310.18515) for details).
- `ppiref_6A_filtered_clustered_04.json` - a non-redundant subset of PPIs that are deduplicated at 0.04 threshold of iDist (a.k.a PPIRef50K).

The latter subset may be the best choice for training machine learning models since it removes biases towards interactions that are overrepresneted in PDB. You can read the subset via the following code (see the corresponding section below for more examples):
```python
from ppiref.split import read_fold
ppi_paths = read_fold('ppiref_6A_filtered_clustered_04', 'whole')
```

Further, this package provides methods to analyze, compare and deduplicate PPIs, as well as to search for similar PPIs in PDB. Please see the examples below.

# PPIRef package

## Installation

If you do not have Python installed, we recommend using the [Miniconda distribution](https://docs.anaconda.com/free/miniconda/). Please install miniconda and create a new environment:
```
conda create -n ppiref python=3.10
conda activate ppiref
```

Then, you can clone and install the complete PPIRef package by running
```
git clone https://github.com/anton-bushuiev/PPIRef.git
cd PPIRef; pip install -e .
```
Alternatively, if you are not interested in using the complete repository, including datasets, external software (`./external`), scripts (`./scripts`) or tests (`./tests`), you can simply install PPIRef as a package without cloning the repository:
```
pip install git+https://github.com/anton-bushuiev/PPIRef.git
```
This is suitable when you are using the PPIRef code as a dependency and will only install Python source code and download basic data files (data splits and example files from `ppiref/data`).

Please see the `external/README.md` directory for the details on how to install the external software for comparing PPIs and calculating buried surface area (BSA).

## Extracting PPIs

The `PPIExtractor` class enables extracting PPIs from PDB files based on inter-atomic distances.

```python
from ppiref.extraction import PPIExtractor
from ppiref.definitions import PPIREF_TEST_DATA_DIR

# Initialize PPI extractor based on 10A contacts between heavy atoms
# Additionally, caluclate buried surface area (BSA) of PPIs (slow)
ppi_dir = PPIREF_TEST_DATA_DIR / 'ppi_dir'
extractor = PPIExtractor(out_dir=ppi_dir, kind='heavy', radius=10., bsa=True)

# Extract all contact-based dimeric PPIs from a PDB file
pdb_file = PPIREF_TEST_DATA_DIR / '1bui.pdb'
extractor.extract(pdb_file)

# Extract a contact-based PPI between two specified chains (dimer)
extractor.extract(pdb_file, partners=['A', 'C'])

# Extract a contact-based PPI between three specified chains (trimer)
extractor.extract(pdb_file, partners=['A', 'B', 'C'])

# Extract a complete dimer complex by setting high expansion radius around interface
ppi_complexes_dir = PPIREF_TEST_DATA_DIR / 'ppi_dir_complexes'
extractor = PPIExtractor(out_dir=ppi_complexes_dir, kind='heavy', radius=6., bsa=False, expansion_radius=1_000_000.)
extractor.extract(pdb_file, partners=['A', 'C'])
```

## Analysing PPIs

After the extraction, one can visualize the PPIs via the PyMOL wrapper and get their statistics.

```python
from ppiref.visualization import PyMOL
pymol = PyMOL()

# Visualize extracted PPIs in PyMOL session + static image
ppi_file = ppi_complexes_dir / 'bu/1bui_A_C.pdb'  # complex
pymol.display_ppi(ppi_file, sticks=False, letters=False, color_by_residues=False)
```

<p align="center">
  <img width="500" src="assets/readme-complex.png"/>
</p>

```python
ppi_file = ppi_dir / 'bu/1bui_A_C.pdb'  # interface
pymol.display_ppi(ppi_file, sticks=True, letters=True, color_by_residues=False)
```

<p align="center">
  <img width="500" src="assets/readme-interface.png"/>
</p>

```python
from ppiref.utils.ppi import PPI

# Get properties of a PPI
ppi = PPI(ppi_file)
ppi.stats
> {'KIND': 'heavy',
>  'EXTRACTION RADIUS': 10.0,
>  'EXPANSION RADIUS': 0.0,
>  'RESOLUTION': 2.65,
>  'STRUCTURE METHOD': 'x-ray diffraction',
>  'DEPOSITION DATE': '1998-09-04',
>  'RELEASE DATE': '1999-09-02',
>  'BSA': 898.378288445}
```

## Comparing PPIs

This package provides wrappers for iAlign and US-align (see `external/README.md`), as well as their scalable approximation iDist (see our [paper](https://arxiv.org/pdf/2310.18515.pdf)) for comparing PPI structures. Additionally it provides a sequence identity comparator to compare PPIs by the corresponding sequences.

```python
from ppiref.comparison import IAlign, USalign, IDist, SequenceIdentityComparator

# Prepare near-duplicate PPIs from Figure 1 in the paper (https://arxiv.org/pdf/2310.18515.pdf)
extractor = PPIExtractor(out_dir=ppi_dir, kind='heavy', radius=6., bsa=False)
extractor.extract(PPIREF_TEST_DATA_DIR / '1p7z.pdb', partners=['A', 'C'])
extractor.extract(PPIREF_TEST_DATA_DIR / '3p9r.pdb', partners=['B', 'D'])
ppis = [ppi_dir / 'p7/1p7z_A_C.pdb', ppi_dir / 'p9/3p9r_B_D.pdb']

# Compare with iAlign (see external/README.md for installation)
# High IS-score and low P-value indicate high similarity
ialign = IAlign()
ialign.compare(*ppis)
> {'PPI0': '1p7z_A_C', 'PPI1': '3p9r_B_D', 'IS-score': 0.95822, 'P-value': 8.22e-67, 'Z-score': 152.167, 'Number of aligned residues': 249, 'Number of aligned contacts': 347, 'RMSD': 0.37, 'Seq identity': 0.992}

# Compare with US-align (see external/README.md for installation)
# High TM-scores in both directions indicate high similarity
usalign = USalign()
usalign.compare(*ppis)
> {'PPI0': '1p7z_A_C', 'PPI1': '3p9r_B_D', 'TM1': 0.992, 'TM2': 0.9965, 'RMSD': 0.3, 'ID1': 0.991, 'ID2': 0.996, 'IDali': 0.998, 'L1': 448, 'L2': 446, 'Lali': 445}

# Compare with iDist
# Low L2 (Euclidean) distance indicates high similarity (below 0.04 is considered near-duplicate for 6A distance interfaces)
idist = IDist()
idist.compare(*ppis)
> {'PPI0': '1p7z_A_C', 'PPI1': '3p9r_B_D', 'L2': 0.0026614179313795114, 'L1': 0.006036636849518753 'Cosine Similarity': 0.999777940667365}

# Compare by maximum pairwise sequence identity
# High sequence identity indicates high similarity
seqid = SequenceIdentityComparator(pdb_dir=PPIREF_TEST_DATA_DIR)  # also requires dir with complete PDB files to extract sequences
seqid.compare(*ppis)
> {'PPI0': '1p7z_A_C', 'PPI1': '3p9r_B_D', 'Maximum pairwise sequence identity': 0.9944979367262724}

# Compare PPIs pairwise with iDist
# (possible with other methods as well but has heavy quadratic complexity)
idist.compare_all_against_all(ppis, ppis)
>            PPI0          PPI1        L2        L1  Cosine Similarity
> 0  1p7z_A_C.pdb  1p7z_A_C.pdb  0.000000  0.000000           1.000000
> 1  1p7z_A_C.pdb  3p9r_B_D.pdb  0.002661  0.006037           0.999778
> 2  3p9r_B_D.pdb  1p7z_A_C.pdb  0.002661  0.006037           0.999778
> 3  3p9r_B_D.pdb  3p9r_B_D.pdb  0.000000  0.000000           1.000000
```

## Deduplicating PPIs

Since iDist is based on a simple vectorization of PPIs, it can be used to deduplicate them based on the Euclidean (L2) distance. The validated threshold used by default is 0.04, suitable for 6A heavy-atom interfaces.

```python
# Embed PPIs with iDist (without comparing)
idist.embed(ppis[0])
idist.embed(ppis[1])
idist.embeddings
> {'3p9r_B_D': array([0.02247098, 0.00224901, 0.03596857, 0.03821186, 0.03148652,
>         0.0404701 , 0.02247443, 0.02697991, 0.01798007, 0.03372832,
>         0.0044929 , 0.0382254 , 0.04607898, 0.01573758, 0.03933556,
>         0.0247265 , 0.0314755 , 0.01123243, 0.00449045, 0.01348696]),
>  '1p7z_A_C': array([0.02237079, 0.00223889, 0.03804226, 0.0391576 , 0.03134589,
>         0.04028928, 0.02237348, 0.02685891, 0.01789972, 0.03357791,
>         0.00447256, 0.0380533 , 0.04587724, 0.01566676, 0.03804791,
>         0.02461776, 0.0313348 , 0.01118212, 0.00447062, 0.01342645])}

# And then query for near-duplicates
# (the query returns the codes of near-duplicate PPIs and their distances to the query embedding)
qeury_emb = idist.embeddings['1p7z_A_C']  # reusing the embedding of the 2nd PPI for example
idist.query(qeury_emb)
> (array([0.        , 0.00346618]),
>  array(['1p7z_A_C', '3p9r_B_D'], dtype=object))

# Or directly deduplicate them
# (the method removes the near-duplicate PPIs from the embeddings)
idist.deduplicate_embeddings()
idist.embeddings
> {'3p9r_B_D': array([0.02247098, 0.00224901, 0.03596857, 0.03821186, 0.03148652,
>         0.0404701 , 0.02247443, 0.02697991, 0.01798007, 0.03372832,
>         0.0044929 , 0.0382254 , 0.04607898, 0.01573758, 0.03933556,
>         0.0247265 , 0.0314755 , 0.01123243, 0.00449045, 0.01348696])}
```

## Finding similar PPIs in Protein Data Bank

The package enables fast search for similar PPIs in PDB based on the interface structure or protein sequence of interest. In the following examples, we will use the same example PPI as in the "Comparing PPIs" section, the KatE homooligomer (1p7z_A_C). Fast search requires precomputed data: iDist embeddings for interface search and MMseqs2 database for sequence search. Thereofore, please download the `ppiref_6A_stats.zip` from [Zenodo](https://zenodo.org/records/10700674), unzip and put it under `PPIRef/ppiref/data/ppiref/ppi_6A_stats` to follow the examples. This can be done automatically by running:
```python
from ppiref.utils.misc import download_from_zenodo
download_from_zenodo('ppi_6A_stats.zip')
```

**By similar interface structure**

One can find PPIs in PDB that are structurally similar to the query PPI. This can be done using the precomputed iDist embeddings. Under the hood, iDist will build an `sklearn` index for all the PPI embeddings and use it to find the neighbors of the query embedding, in the near-duplicate radius (0.04 by default, which is validated for 6A interfaces).

```python
from ppiref.definitions import PPIREF_DATA_DIR, PPIREF_TEST_DATA_DIR
from ppiref.comparison import IDist

# Initialize IDist and read embeddings for all PPI interfaces in PPIRef (i.e., all PPIs in PDB)
idist = IDist()
idist.read_embeddings(PPIREF_DATA_DIR / 'ppiref/ppi_6A_stats/idist_emb.csv', dropna=True)

# Embed your query PPI (without necessarily storing in the iDist embeddings database)
ppi_dir = PPIREF_TEST_DATA_DIR / 'ppi_dir'
query_ppi_path = ppi_dir / 'p7/1p7z_A_C.pdb'
query_embedding = idist.embed(query_ppi_path, store=False)

# Query for similar PPIs
dists, ppi_ids = idist.query(query_embedding)

# Print the 10 most similar PPIs sorted by iDist distances
for dist, ppi_id in list(zip(dists, ppi_ids))[:10]:
    print(f'{ppi_id} {dist:.4f}')
> 1p7z_A_C 0.0000
> 4env_A_C 0.0017
> 1p80_A_C 0.0024
> 1p7y_A_C 0.0024
> 4enu_A_C 0.0030
> 4ens_A_C 0.0030
> 3p9s_B_D 0.0034
> 1p7z_B_D 0.0035
> 4ent_A_C 0.0035
> 3p9r_B_D 0.0035
```

The top hit is the same interaction as the query one, 1p7z_A_C, with the zero distance. The other hits are near-duplicates of the query PPI (compare, for example, the [1p7z_A_C](https://www.rcsb.org/structure/1P7Z) and [4env_A_C](https://www.rcsb.org/structure/4env) entries).

**By similar sequeunce**

One can also find PPIs in PDB that involve sequences similar to the one of interest. This can be done using the prepared [MMseqs2](https://github.com/soedinglab/mmseqs2) database. Install MMseqs2 according to the official documentation and then you can use the wrapper as below. Under the hood, the wrapper will use the `mmseqs2 easy-search` with default parameters.

```python
from ppiref.comparison import MMSeqs2PPIRetriever

# Initialize the wrapper for MMseqs2 database to store all sequences from PPIRef
mmseqs2 = MMSeqs2PPIRetriever(PPIREF_DATA_DIR / 'ppiref/ppi_6A_stats/mmseqs_db/db')

# Prepare your fasta file (for example by downloading from Uniprot)
query_path = PPIREF_TEST_DATA_DIR / '1p7z_A.fasta'

# Query the MMseqs2 database for PPIs involving a sequence similar to the query sequence
# (returns triples (PPI id, sequence similarity, partner similar to query sequence))
seq_sims, ppi_ids, partners =  mmseqs2.query(query_path)
for ppi_id, seq_sim, partner in list(zip(ppi_ids, seq_sims, partners))[:10]:
    print(f'{ppi_id} {seq_sim} {partner}')
> 4bfl_A_B 1.0 B
> 4bfl_B_C 1.0 B
> 4bfl_B_D 1.0 B
> 1iph_C_D 1.0 C
> 1iph_B_C 1.0 C
> 1iph_A_C 1.0 C
> 1gge_A_C 1.0 A
> 1gge_A_B 1.0 A
> 1gge_A_D 1.0 A
> 3vu3_A_H 1.0 A
```

The retrieved PPIs involve the same sequence (100% identity) as the query one (chain A from 1p7z), and 
they are proper PPIs, as they satisfy the PPIRef filtering criteria (see "PPIRef dataset" section). One can also validate that the retrieved PPIs include the originating PPI, as well as, for example, the top iDist hit from above:

```python
'1p7z_A_C' in ppi_ids and '4env_A_C' in ppi_ids
> True
```

## Splitting and subsetting PPIs

The package provides a unified approach to storing and processing data splits and other subsets of PPIs.

```python
from ppiref.split import read_split, read_fold, read_split_source, write_split

# Read PPI codes for all 6A interfaces from PPIRef
# (from a prepared .json file in ./ppiref/data/splits)
split = read_split('ppiref_6A_raw', full_paths=False)
split['whole'][:3]
> ['6q2a_F_O', '3q2k_I_P', '8q2k_A_C']

# Read location of .pdb files of the split
read_split_source('ppiref_6A_raw')
> PosixPath('/scratch/project/open-26-23/antonb/PPIRef/ppiref/data/ppiref/ppi_6A')

# Or direcrly read the fold with full paths (whole dataset in this case)
fold = read_fold('ppiref_6A_raw', 'whole', full_paths=True)
fold[0]
> PPIPath('/scratch/project/open-26-23/antonb/tmp/PPIRef/ppiref/data/ppiref/ppi_6A/q2/6q2a_F_O.pdb')
```

Examples of predefined splits:

```python
# PPIRef50K used to train PPIformer
fold = read_fold('ppiref_10A_filtered_clustered_03', 'whole', full_paths=False)
fold[:3]
> ['4q2p_A_B', '3q2s_A_D', '6q2v_B_E']

# Test set from non-leaking SKEMPI v2.0
fold = read_fold('skempi2_iclr24_split', 'test', full_paths=False)
fold[:3]
> ['1B3S_A_D', '1B2U_A_D', '1BRS_A_D']

# DIPS set used to train and validate EquiDock and DiffDock-PP
fold = read_fold('dips_equidock', 'train+val', full_paths=False)
> ['1v6j_A_D', '2v6a_A_L', '2v6a_B_O']

# Write your own split (this will also run simple sanity checks)
split = {'train': ['1p7z_A_C'], 'test': ['3p9r_B_D', '1p7z_A_C']}
write_split('demo_split', source=PPIREF_TEST_DATA_DIR / 'ppi_dir', folds=split)
> UserWarning: Folds train and test are not disjoint.
> UserWarning: Split is not complete: 2 of 5 PPIs contained.
```

# TODO

Technical

- [x] PPIRef (6A interfaces) on Zenodo
- [ ] PPIRef (10A interfaces) on Zenodo (expected in June 2024)
- [ ] PPIRef version based on the [PDB-REDO database](https://pdb-redo.eu/) for higher-quality side chains in the structures (expected in June 2024)
- [ ] Docstrings

Enhancements

- [ ] Cluster all PPIs to sample from clusters rather than removing near duplicates completely (similar to UniRef seeds)
- [ ] Add RASA values to classify residues according to [Levy 2010](https://pubmed.ncbi.nlm.nih.gov/20868694/)
- [ ] Classify PPIs according to [Ofran2003](https://pubmed.ncbi.nlm.nih.gov/12488102/)

# References

If you find this repository useful, please cite our paper or the corresponding external software (see `external/README.md`).
```
@article{
  bushuiev2024learning,
  title={Learning to design protein-protein interactions with enhanced generalization},
  author={Anton Bushuiev and Roman Bushuiev and Petr Kouba and Anatolii Filkin and Marketa Gabrielova and Michal Gabriel and Jiri Sedlar and Tomas Pluskal and Jiri Damborsky and Stanislav Mazurenko and Josef Sivic},
  booktitle={The Twelfth International Conference on Learning Representations},
  year={2024}
}
```
