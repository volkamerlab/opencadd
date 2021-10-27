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

<!---
A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.
-->

Protein kinases are involved in most aspects of cell life due to their role in signal transduction. Dysregulated kinases can cause severe diseases such as cancer, inflammatory and neurodegenerative diseases, which has made them a frequent target in drug discovery for the last decades [@Cohen:2001].
The immense research on kinases has led to an increasing amount of kinase resources [@Kooistra:2017].
Among them is the KLIFS database, which focuses on storing and analyzing structural data on the binding of drugs and other small molecules to kinases [@Kanev:2021].
Convenient integration of the KLIFS data into workflows is the aim of the OpenCADD-KLIFS Python module presented here.

# Statement of need

<!---
A Statement of Need section that clearly illustrates the research purpose of the software.
-->

OpenCADD is a collection of Python modules for structural cheminformatics; OpenCADD-KLIFS is one of them (``opencadd.databases.klifs``).
This module offers access to KLIFS data such as information about kinases, structures, ligands, 
interaction fingerprints, and bioactivities. 
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or remotely from the KLIFS webserver. 
The tool provides identical APIs for the remote and local queries and streamlines all output into 
standardized Pandas [@pandas] DataFrames to allow for an easy and quick downstream data manipulation (Figure \autoref{fig:opencadd_klifs_toc}). This Pandas-focused setup is ideal to work with Jupyter notebooks [@jupyterhub; @Kluyver:2016]. 


![OpenCADD-KLIFS fetches KLIFS [@Kanev:2021] data offline from a KLIFS download or online from the KLIFS database and formats the output in user-friendly Pandas [@pandas] DataFrames.\label{fig:opencadd_klifs_toc}](opencadd_klifs_toc.png)

The KLIFS database offers a REST API including an OpenAPI specification. OpenCADD-KLIFS uses bravado [@bravado] to dynamically generate a Python client based on the OpenAPI definitions and add wrappers to enable the following functionalities:

- A session is set up, which allows access to various KLIFS data sources by different identifiers with the API ``session.data_source.by_identifier``; for example ``session.structures.by_kinase_name`` fetches information on all structures for a query kinase. Data sources currently include kinases, ligands, structures, drugs, pockets, bioactivities, structural conformations, modified residues, and coordinates; identifiers refer to kinase names, PDB IDs, KLIFS IDs, and more.
- Query results obtained from the remote KLIFS webserver are streamlined with those obtained from a local KLIFS download using the same API.
- All results obtained with bravado are formatted as Pandas DataFrames with standardized column names, data types, and handling of missing data.
- Structural files deposited on KLIFS include full complexes or selections such as proteins, pockets, ligands, and more. These files can be downloaded to disc or loaded via biopandas [@biopandas] or the RDKit [@rdkit].

OpenCADD-KLIFS is especially useful whenever users are interested in multiple or more complex queries such as "fetch all interaction profiles of kinases bound to the drug Gefitinib". Formatting the output as DataFrames facilitates subsequent filtering steps and DataFrame merges in case multiple KLIFS datasets need to be combined.
OpenCADD-KLIFS is currently used in several projects from the Volkamer Lab [@volkamerlab] including TeachOpenCADD [@teachopencadd], OpenCADD-pocket [@opencaddpocket], KiSSim [@kissim], KinoML [@kinoml], and PLIPify [@plifify].

# Acknowledgements

We thank Albert Kooistra for his help with questions and wishes regarding the KLIFS database, and David Schaller for his feedback on the OpenCADD-KLIFS module.
We acknowledge the contributors involved in software programs and packages used by OpenCADD-KLIFS, such as bravado, RDKit, Pandas, Jupyter, and Pytest, and Sphinx. 

# References