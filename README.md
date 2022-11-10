# CellPacking

This is a flat monolayer simulation package for the article: 
### Curvature-induced cell rearrangements in biological tissues

Yuting Lou<sup>1</sup>, Jean-Francois Rupprecht<sup>1,2</sup>, Sophie Theis<sup>3</sup>, Tetsuya Hiraiwa<sup>1</sup>, and Timothy E Saunders<sup>1,3</sup>

<sup>1</sup>Mechanobiology Institute, National University of Singapore  
<sup>2</sup>Aix Marseille Université, Université de Toulon, CNRS,  
Centre de Physique Théorique, Turing Centre for Living Systems, Marseille, France  
<sup>3</sup>Warwick Medical School, University of Warwick, Coventry, United Kingdom


## Installation
This package is based on the [`tyssus`](https://tyssue.readthedocs.org) library and its dependencies. 

The recommanded installation route is to use the `conda` package manager. You can get a `conda` distribution for your OS at https://www.anaconda.com/download . Make sure to choose a python 3.6 version. Once you have installed conda, you can install tyssue with:

```bash
$ conda install -c conda-forge tyssue
```

You can then download and install CellPacking from github:

- with git:

```bash
$ git clone https://github.com/TimSaundersLab/CellPacking.git
$ cd invagination
$ python setup.py install
```

- or by downloading https://github.com/TimSaundersLab/CellPacking/archive/master.zip ,  uncompressing the archive and running `python setup.py install` in the root directory.

## Licence

This work is free software, published under the MPLv2 licence, see LICENCE for details.


&copy; The article authors -- all rights reserved
