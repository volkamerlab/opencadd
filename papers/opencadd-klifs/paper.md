---
title: 'OpenCADD-KLIFS: A Python package to fetch kinase data from the KLIFS database'
tags:
  - Python
  - KLIFS
  - kinase
authors:
  - name: Dominique Sydow^[corresponding author]
    orcid: 0000-0003-4205-8705
    affiliation: 1
  - name: Jaime Rodríguez-Guerra
    affiliation: 2
  - name: Andrea Volkamer
    affiliation: 1
affiliations:
 - name: _In Silico_ Toxicology and Structural Bioinformatics, Institute of Physiology, Charité – Universitätsmedizin Berlin, corporate member of Freie Universität Berlin and Humboldt-Universität zu Berlin, Augustenburger Platz 1, 13353 Berlin, Germany
   index: 1
 - name: Institution Name
   index: 2
date: 27 October 2021
bibliography: paper.bib
---

# Summary

Protein kinases are involved in most aspects of cell life due to their role in signal transduction. Dysregulated kinases can cause severe diseases such as cancer, inflammatory and neurodegenerative diseases, which has made them a frequent target in drug discovery for the last decades [@Cohen:2001].
The immense research on kinases has led to an increasing amount of kinase resources [@Kooistra:2017].
Among them is the KLIFS database, which focuses on storing and analyzing structural data on kinases and interacting drugs and other small molecules [@Kanev:2021].
The OpenCADD-KLIFS Python module offers a convenient integration of the KLIFS data into workflows to facilitate computational kinase research.

# Statement of need

[OpenCADD-KLIFS](https://opencadd.readthedocs.io/en/latest/databases_klifs.html) (``opencadd.databases.klifs``) is a part of the [OpenCADD](https://opencadd.readthedocs.io/) package, a collection of Python modules for structural cheminformatics.
This module offers access to KLIFS data [@Kanev:2021; @klifs_website] such as information about kinases, structures, ligands, 
interaction fingerprints, and bioactivities. 
KLIFS thereby focuses especially on the ATP binding site, defined as a set of 85 residues and aligned across all structures using a multiple sequence alignment (MSA) [vanLinden:2014].
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or remotely from the KLIFS webserver. 
The presented module provides identical APIs for the remote and local queries for KLIFS data and streamlines all output into 
standardized Pandas DataFrames [@pandas] to allow for easy and quick downstream data analyses (Figure \autoref{fig:opencadd_klifs_toc}). This Pandas-focused setup is ideal to work with in Jupyter notebooks [@jupyterhub; @Kluyver:2016]. 


![OpenCADD-KLIFS fetches KLIFS data [@Kanev:2021] offline from a KLIFS download or online from the KLIFS database and formats the output in user-friendly Pandas DataFrames [@pandas].\label{fig:opencadd_klifs_toc}](opencadd_klifs_toc.png)

The KLIFS database offers a REST API including an OpenAPI specification [@klifs_swagger]. Our module OpenCADD-KLIFS uses bravado [@bravado] to dynamically generate a Python client based on these OpenAPI definitions and adds wrappers to enable the following functionalities:

- A session is set up, which allows access to various KLIFS *data sources* by different *identifiers* with the API ``session.data_source.by_identifier``. *Data sources* currently include kinases, structures and annotated conformations, modified residues, pockets, ligands, drugs, and bioactivities; *identifiers* refer to kinase names, PDB IDs, KLIFS IDs, and more.
For example, ``session.structures.by_kinase_name`` fetches information on all structures for a query kinase.
- Query results obtained from the remote KLIFS webserver are streamlined with those obtained from a local KLIFS download using the same API.
- All results obtained with bravado are formatted as Pandas DataFrames with standardized column names, data types, and handling of missing data.
- Files with the structural 3D coordinates deposited on KLIFS include full complexes or selections such as proteins, pockets, ligands, and more. These files can be downloaded to disc or loaded via biopandas [@biopandas] or RDKit [@rdkit].

OpenCADD-KLIFS is especially convenient whenever users are interested in multiple or more complex queries such as "fetching all structures for the kinase EGFR in the DFG-in conformation" or "fetching the measured bioactivity profiles for all ligands that are structurally resolved in complex with EGFR". Formatting the output as DataFrames facilitates subsequent filtering steps and DataFrame merges in case multiple KLIFS datasets need to be combined.
OpenCADD-KLIFS is currently used in several projects from the Volkamer Lab [@volkamerlab] including TeachOpenCADD [@teachopencadd], OpenCADD-pocket [@opencadd_pocket], KiSSim [@kissim], KinoML [@kinoml], and PLIPify [@plipify].
For example, OpenCADD-KLIFS is applied in a [TeachOpenCADD tutorial](https://projects.volkamerlab.org/teachopencadd/talktorials/T012_query_klifs.html) to demonstrate how to fetch all kinase-ligand interaction profiles for all available EGFR kinase structures to visualize the per-residue interaction types and frequencies with only a few lines of code.

# Acknowledgements

We thank the whole KLIFS team for providing such a great kinase resource with an easy-to-use API and especially Albert Kooistra for his help with questions and wishes regarding the KLIFS database. 
We thank David Schaller for his feedback on the OpenCADD-KLIFS module.
We acknowledge the contributors involved in software programs and packages used by OpenCADD-KLIFS, such as bravado, RDKit, Pandas, Jupyter, and Pytest, and Sphinx. 

# References