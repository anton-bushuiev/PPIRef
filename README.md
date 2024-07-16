<!-- <div align="center"> -->

# PPIRef

[![Documentation badge](https://img.shields.io/badge/docs-latest-brightgreen.svg)](https://ppiref.readthedocs.io/en/latest/?badge=latest)
[![arXiv badge](https://img.shields.io/badge/arXiv-2310.18515-b31b1b.svg)](https://arxiv.org/abs/2310.18515)
[![Zenodo badge](https://zenodo.org/badge/DOI/10.5281/zenodo.12746715.svg)](https://doi.org/10.5281/zenodo.12746715)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python package](https://github.com/anton-bushuiev/PPIRef/actions/workflows/python-package.yml/badge.svg)](https://github.com/anton-bushuiev/PPIRef/actions/workflows/python-package.yml)
![Python Versions](https://img.shields.io/badge/Python-3.9%20%7C%203.10%20%7C%203.11-green.svg)

<!-- </div> -->

<p align="center">
  <img src="https://raw.githubusercontent.com/anton-bushuiev/PPIRef/f967861bd665e36d13dec493f054a1b2a9dd5538/assets/readme-dimer-close-up.png"/>
</p>

[PPIRef](https://github.com/anton-bushuiev/PPIRef/tree/main) is a Python package for working with 3D structures of protein-protein interactions (PPIs). It is based on the PPIRef dataset, comprising all PPIs from the Protein Data Bank (PDB). The package aims to provide standard data and tools for machine learning and data science applications involving protein-protein interaction structures. PPIRef includes the following functionalities:

- ⭐ **Extracting** protein-protein interfaces from .pdb files.
- ⭐ **Visualizing and analyzing** the properties of PPIs.
- ⭐ **Comparing, deduplicating and clustering** PPI interfaces.
- ⭐ **Retrieving** similar PPIs from PDB by similar interface structure or sequence.
- ⭐ **Downloading, splitting and subsampling** prepared PPIs for machine learning applications.

Please see the [documentation](https://ppiref.readthedocs.io/en/latest/) for usage examples and API reference. See also [our paper](https://arxiv.org/abs/2310.18515) for additional details.

## Quick start 🚀

Install the PPIRef package.

```bash
conda create -n ppiref python=3.10
conda activate ppiref
git clone https://github.com/anton-bushuiev/PPIRef.git
cd PPIRef; pip install -e .
```

Download the dataset using the package (in Python).

```python
from ppiref.utils.misc import download_from_zenodo
from ppiref.split import read_fold
from ppiref.utils.ppi import PPI
download_from_zenodo('ppi_6A.zip')  # or for example 'pdb_redo_ppi_10A.zip' for all 10-Angstrom PPIs from PDB-REDO
> Downloading: 100%|██████████| 6.94G/6.94G [10:19<00:00, 11.2MiB/s]
> Extracting: 100%|██████████| 831382/831382 [02:36<00:00, 5313.49files/s]
```

Read the data fold/subset you need (whole PPIRef50K in the example).

```python
ppi_paths = read_fold('ppiref_6A_filtered_clustered_04', 'whole')
print('Dataset size:', len(ppi_paths))
> Dataset size: 51755
```

Now you are ready to work with the PPIRef dataset! Example of a sample:

```python
ppi = PPI(ppi_paths[0])
print('Path:', ppi.path)
print('Statistics:', ppi.stats)
ppi.visualize()
> Path: /Users/anton/dev/PPIRef/ppiref/data/ppiref/ppi_6A/hc/3hch_A_B.pdb
> Statistics: 
> {'KIND': 'heavy',
>  'EXTRACTION RADIUS': 6.0,
>  'EXPANSION RADIUS': 0.0,
>  'RESOLUTION': 2.1,
>  'STRUCTURE METHOD': 'x-ray diffraction',
>  'DEPOSITION DATE': '2009-05-06',
>  'RELEASE DATE': '2009-10-13',
>  'BSA': 682.5337386399999}
```

<p align="center">
  <img width=500px src="https://raw.githubusercontent.com/anton-bushuiev/PPIRef/5fca49ecd0e776a362e6f8dc090f14432b6efbd6/assets/3hch_A_B.png"/>
</p>

Further, the PPIRef package provides utilities for comparing, deduplicating, and clustering PPI interfaces, as well as for retrieving similar PPIs from PDB by similar interface structure or sequence. Please see the [documentation](https://ppiref.readthedocs.io/en/latest/) for more details.

## TODO

The repository is under development. Please do not hesitate to contact us or create an issue/PR if you have any questions or suggestions ✌️.

**Technical**

- [x] PPIRef (6A interfaces) on Zenodo
- [x] PPIRef (10A interfaces) on Zenodo (expected in June 2024)
- [x] PPIRef version based on the [PDB-REDO database](https://pdb-redo.eu/) for higher-quality side chains in the structures (expected in June 2024)
- [x] Docstrings

**Enhancements**

- [ ] Cluster all PPIs to sample from clusters rather than removing near duplicates completely (similar to UniRef seeds)
- [ ] Add RASA values to classify residues according to [Levy 2010](https://pubmed.ncbi.nlm.nih.gov/20868694/)
- [ ] Classify PPIs according to [Ofran2003](https://pubmed.ncbi.nlm.nih.gov/12488102/)

## References

If you find this repository useful, please cite our paper:

```bibtex
@article{bushuiev2024learning,
  title={Learning to design protein-protein interactions with enhanced generalization},
  author={Anton Bushuiev and Roman Bushuiev and Petr Kouba and Anatolii Filkin and Marketa Gabrielova and Michal Gabriel and Jiri Sedlar and Tomas Pluskal and Jiri Damborsky and Stanislav Mazurenko and Josef Sivic},
  booktitle={ICLR 2024 (The Twelfth International Conference on Learning Representations)},
  url={https://doi.org/10.48550/arXiv.2310.18515},
  year={2024}
}
```

If relevant, please also cite the corresponding paper on data leakage in protein interaction benchmarks:

```bibtex
@article{bushuiev2024revealing,
  title={Revealing data leakage in protein interaction benchmarks},
  author={Anton Bushuiev and Roman Bushuiev and Jiri Sedlar and Tomas Pluskal and Jiri Damborsky and Stanislav Mazurenko and Josef Sivic},
  booktitle={ICLR 2024 Workshop on Generative and Experimental Perspectives for Biomolecular Design},
  url={https://doi.org/10.48550/arXiv.2404.10457},
  year={2024}
}
```

If you find any of the external software useful, please cite the corresponding papers (see `PPIRef/external/README.md`).
