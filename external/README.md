# External software

## iAlign
The `ppiref.comparison.IAlign` class is a Python wrapper for the [iAlign software](https://sites.gatech.edu/cssb/ialign/). iAlign is an alignment-based algorithm for the structural comparison of protein-protein interfaces. Essentially, iAlign is an adaptation of TM-score to the domain of protein-protein interactions.

To use the wrapper, please download the official Perl source code from the [iAlign website](https://sites.gatech.edu/cssb/ialign/) and place it in this directory. Alternatively, you can place iAlign under a different location but change the `ppiref.defitions.IALIGN_PATH`. After installation, the resulting directory structure may look like this:

```
iAlign
├── README.md
├── bin
└── example

3 directories, 1 file
```

If you find iAlign useful, please cite the original paper:
```
@article{gao2010ialign,
  title={iAlign: a method for the structural comparison of protein--protein interfaces},
  author={Gao, Mu and Skolnick, Jeffrey},
  journal={Bioinformatics},
  volume={26},
  number={18},
  pages={2259--2265},
  year={2010},
  publisher={Oxford University Press}
}
```

## US-align

The `ppiref.comparison.USAlign` class is a Python wrapper for the [US-align software](https://zhanggroup.org/US-align/#:~:text=US%2Dalign%20standalone%20program%20download) to compare PPI structures. US-align is a more recent adaption of TM-score and is origianlly designed for the universal comparison of different kinds of macromolecules.

To use, the wrapper, please download the official compiled C++ executable from the [US-align website](https://zhanggroup.org/US-align/#:~:text=US%2Dalign%20standalone%20program%20download) and place it in the `ppiref/external` directory. Alternatively, you can place US-align under a different location but change the `ppiref.defitions.USALIGN_PATH`. The resulting directory structure may look like this:
```
USalign
└── USalign

1 directories, 1 files
```

If you find US-align useful, please cite the original paper:
```
@article{zhang2022us,
  title={US-align: universal structure alignments of proteins, nucleic acids, and macromolecular complexes},
  author={Zhang, Chengxin and Shine, Morgan and Pyle, Anna Marie and Zhang, Yang},
  journal={Nature methods},
  volume={19},
  number={9},
  pages={1109--1115},
  year={2022},
  publisher={Nature Publishing Group US New York}
}
```

## dr_sasa
The `ppiref.surface.DR_SASA` class is a Python wrapper for the [dr_sasa software](https://github.com/nioroso-x3/dr_sasa_n) for calculating buried surface area (BSA) of a PPI interface. 

In order to use it, build the C++ source code according to the [original doctumentation](https://github.com/nioroso-x3/dr_sasa_n#compiling). The resulting executable should match the `ppiref.defitions.DR_SASA_PATH`. By default the path expects building dr_sasa in this (`PPIRef/external`) directory. The resulting directory structure may look like this:
```
dr_sasa_n
├── CMakeLists.txt
├── INSTALL
├── LICENSE
├── README.md
├── build
├── doc
├── examples
├── src
└── utils

6 directories, 4 files
```

If you find dr_sasa useful, please cite the original paper:
```
@article{ribeiro2019calculation,
  title={Calculation of accurate interatomic contact surface areas for the quantitative analysis of non-bonded molecular interactions},
  author={Ribeiro, Judemir and R{\'\i}os-Vera, Carlos and Melo, Francisco and Sch{\"u}ller, Andreas},
  journal={Bioinformatics},
  volume={35},
  number={18},
  pages={3499--3501},
  year={2019},
  publisher={Oxford University Press}
}
```
