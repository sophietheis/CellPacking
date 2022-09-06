# CellPacking

This is a flat monolayer simulation package for the article: 
### Curvature-induced cell rearrangements in biological tissues

Yuting Lou, Jean-Francois Rupprecht, Tetsuya Hiraiwa, and Timothy E Saunders


Try it with my binder by clicking the badge bellow:
[![Binder](https://mybinder.org/badge_logo.svg)](PUT LINK HERE OF A NOTEBOOK)

## Installation
This package is based on the [`tyssus`](https://tyssue.readthedocs.org) library and its dependencies. 

The recommanded installation route is to use the `conda` package manager. You can get a `conda` distribution for your OS at https://www.anaconda.com/download . Make sure to choose a python 3.6 version. Once you have installed conda, you can install tyssue with:

```bash
$ conda install -c conda-forge tyssue
```

You can then download and install CellPacking from github:

- with git:

```bash
$ git clone https://github.com/sophietheis/CellPacking.git
$ cd invagination
$ python setup.py install
```

- or by downloading https://github.com/sophietheis/CellPacking/archive/master.zip ,  uncompressing the archive and running `python setup.py install` in the root directory.

## Licence

This work is free software, published under the MPLv2 licence, see LICENCE for details.


&copy; The article authors -- all rights reserved
