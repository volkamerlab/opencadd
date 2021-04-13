# OpenCADD

[//]: # "Badges"

[![GH Actions Status](https://github.com/volkamerlab/opencadd/workflows/CI/badge.svg)](https://github.com/volkamerlab/opencadd/actions?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/volkamerlab/opencadd/branch/master/graph/badge.svg)](https://codecov.io/gh/volkamerlab/opencadd/branch/master)
[![Documentation Status](https://readthedocs.org/projects/opencadd/badge/?version=latest)](https://opencadd.readthedocs.io)

![OpenCADD](/docs/_static/opencadd.png)

A Python library for structural cheminformatics.

## Overview

> This library is still in early stages of development.

- `compounds.standardization`: standardize chemical records.
- `databases.klifs`: utilities to query the KLIFS database, offline or online.
- `io`: read and write molecules from/to files.
- `structure.pocket`: identification and analysis of protein (sub)pockets.
- `structure.superposition` (formerly `superposer`): superimpose macromolecules using sequence and structural information.

## Documentation

The Documentation will be available soon.

## License

`opencadd` is free software and is licensed under the MIT license. Copyright (c) 2020, Volkamer Lab

## Authors

`opencadd` is the cumulative work of several members of the [Volkamer Lab](https://volkamerlab.org), as well as contributions from students that have participated in our lab. In no particular order:

- Jaime Rodríguez-Guerra, PhD
- Dominique Sydow
- Dennis Köser, Annie Pham, Enes Kurnaz, Julian Pipart (structural superposition, 2020)
- Allen Dumler (standardizer, 2021)

# Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
