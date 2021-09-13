.. opencadd documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OpenCADD's documentation!
=========================================================

OpenCADD is a Python package for structural cheminformatics!

.. image::
   https://github.com/volkamerlab/opencadd/workflows/CI/badge.svg
   :target: https://github.com/volkamerlab/opencadd/actions
.. image::
   https://codecov.io/gh/volkamerlab/opencadd/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/volkamerlab/opencadd/branch/master
.. image::
   https://readthedocs.org/projects/opencadd/badge/?version=latest
   :target: https://opencadd.readthedocs.io/en/latest/
.. image::
   https://img.shields.io/conda/vn/conda-forge/opencadd.svg
   :target: https://anaconda.org/conda-forge/opencadd
.. image::
   https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT

.. raw:: html

   <p align="center">
   <img src="_static/opencadd.png" alt="Subpocket-based structural fingerprint for kinase pockets" width="600"/>
   </p>

.. toctree::
   :maxdepth: 1
   :caption: User guide

   installing

.. toctree::
   :maxdepth: 1
   :caption: Input/output formats

   io
   tutorials/io

.. toctree::
   :maxdepth: 1
   :caption: Structure: Superposition

   superposition
   tutorials/mda
   tutorials/mmligner
   tutorials/theseus

.. toctree::
   :maxdepth: 1
   :caption: Structure: Pocket

   structure_pocket
   tutorials/structure_pocket

.. toctree::
   :maxdepth: 1
   :caption: Databases: KLIFS

   databases_klifs
   tutorials/databases_klifs

.. toctree::
   :maxdepth: 1
   :caption: Developers

   developers
   api




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
