# OpenCADD

[//]: # "Badges"

[![GH Actions Status](https://github.com/volkamerlab/opencadd/workflows/CI/badge.svg)](https://github.com/volkamerlab/opencadd/actions?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/volkamerlab/opencadd/branch/master/graph/badge.svg)](https://codecov.io/gh/volkamerlab/opencadd/branch/master)
[![Documentation Status](https://readthedocs.org/projects/opencadd/badge/?version=latest)](https://opencadd.readthedocs.io)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/opencadd.svg)](https://anaconda.org/conda-forge/opencadd)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/237947037.svg)](https://zenodo.org/badge/latestdoi/237947037)

![OpenCADD](/docs/_static/opencadd.png)

A Python library for structural cheminformatics.

## Overview

> Some modules of this library are still in early stages of development as indicated below.

- `databases.klifs`: utilities to query the KLIFS database, offline or online.
- :construction: `io`: read and write molecules from/to files.
- :construction: `structure.pocket`: identification and analysis of protein (sub)pockets.
- :construction: `structure.superposition` (formerly `superposer`): superimpose macromolecules using sequence and structural information.

## Documentation

The documentation is available [here](https://opencadd.readthedocs.io/en/latest/).

## License

`opencadd` is free software and is licensed under the MIT license. Copyright (c) 2020, Volkamer Lab

## Authors

`opencadd` is the cumulative work of several members of the [Volkamer Lab](https://volkamerlab.org), as well as contributions from students that have participated in our lab. In no particular order:

- Jaime Rodríguez-Guerra, PhD
- Dominique Sydow
- Dennis Köser, Annie Pham, Enes Kurnaz, Julian Pipart (structural superposition, 2020)

# Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
