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
The immense research on kinases has lead to an increasing amount of kinase resources [@Kooistra:2017].
Among them is the KLIFS database, which focuses on storing and analyzing structural data on the binding of drugs and other small molecules to kinases [@Kanev:2021].
Convienent integration of the KLIFS data into workflows is the aim of the OpenCADD-KLIFS Python module presented here.

# Statement of need

<!---
A Statement of Need section that clearly illustrates the research purpose of the software.
-->

OpenCADD is a collection of Python modules for structural cheminformatics; OpenCADD-KLIFS is one of them (``opencadd.databases.klifs``).
This module offers access to KLIFS data such as information about kinases, structures, ligands, 
interaction fingerprints, and bioactivities. 
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or remotely from the KLIFS webserver. 
The tool provides identical APIs for the remote and local queries and streamlines all output into 
standardized Pandas DataFrames to allow for an easy and quick downstream data manipulation (Figure \autoref{fig:opencadd_klifs_toc}).


![OpenCADD-KLIFS fetches KLIFS [@Kanev:2021] data offline from a KLIFS download or online from the KLIFS database and formats the output in user-friendly Pandas [@pandas] DataFrames.\label{fig:opencadd_klifs_toc}](opencadd_klifs_toc.png)

The KLIFS database offers a REST API including an OpenAPI specification. We use bravado [@bravado] to dynamically generate a Python client based on the OpenAPI definitions. With OpenCADD-KLIFS, we offer a Python wrapper around this client to add the following functionalities:
- Format the query results obtained with bravado into Pandas [@pandas] DataFrames with standardized column names and data types. 
- Streamline the query results obtained from the KLIFS webserver with those obtained from a local KLIFS download.
-  TBA

OpenCADD-KLIFS is especially useful whenever users are interested in multiple or more complex queries (e.g. fetch all interaction profiles of kinases bound to Gefitinib). The output Pandas DataFrames allow easy filtering steps and DataFrame merges in case of mulitple queries for different data sources.
OpenCADD-KLIFS is currently used in several projects from the Volkamer Lab [@volkamerlab] including TeachOpenCADD [@teachopencadd], OpenCADD-pocket [@opencaddpocket], and KiSSim [@kissim].

TODO
- Comment on structure file download
- Great for Jupyter notebooks

# Acknowledgements

We thank Albert Kooistra for his help with questions and wishes regarding the KLIFS database.
We acknowledge the contributors involved in software programs and packages used by OpenCADD, such as bravado, RDKit, Pandas, Jupyter, and Pytest, and Sphinx. 

# References