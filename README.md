# StructuralAlignment

[//]: # (Badges)
[![GH Actions Status](https://github.com/volkamerlab/structuralalignment/workflows/CI/badge.svg)](https://github.com/volkamerlab/structuralalignment/actions?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/volkamerlab/StructuralAlignment/branch/master/graph/badge.svg)](https://codecov.io/gh/volkamerlab/StructuralAlignment/branch/master)
[![Documentation Status](https://readthedocs.org/projects/structural-alignment/badge/?version=latest)](https://structural-alignment.readthedocs.io/en/latest/?badge=latest)

A Python library for molecular structural alignment and superposition

## Documentation

The Documentation is available at [readthedocs](https://structural-alignment.readthedocs.io/en/latest/ "Read the Docs").

## Features

The structural alignment and superposition can be performed with different methods like:

* mmligner
* theseus
* Chimera-Matchmaker

## How to install the library

coming soon
<!-- `conda install ...` -->

## Requirements

* Python 3.7 or later

## How to use the library

With this library you can superpose 2 or more proteinstructures.\
If the input is more than 2 structures, the first will be considered the target structure and every other structure will be superposed with this one thereby creating multiple structural alignments.
The sequence length has to be the same.
More detailed instructions can be found in the [documentation](https://structural-alignment.readthedocs.io/en/latest/ "documentation").

<!-- need to add how to use it -->

## Examples

`structuralalignment 1lol 5XME --method mmligner --method-options something`
<!-- need to add examples -->

## Copyright

Copyright (c) 2020, Volkamer Lab

## License

`structuralalignment` is free software and is licensed under the MIT license.

## Authors

* Jaime Rodríguez-Guerra <jaime.rodriguez@charite.de>
* Dennis Köser
* Annie Pham
* Enes Kuni
* Julian Pipart

## Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
