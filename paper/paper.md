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
    orcid: 0000-0001-8974-1566
    affiliation: 1
  - name: Andrea Volkamer
    affiliation: 1
    orcid: 0000-0002-3760-580X
affiliations:
 - name: _In Silico_ Toxicology and Structural Bioinformatics, Institute of Physiology, Charité – Universitätsmedizin Berlin, corporate member of Freie Universität Berlin and Humboldt-Universität zu Berlin, Augustenburger Platz 1, 13353 Berlin, Germany
   index: 1
date: 27 October 2021
bibliography: paper.bib
---

# Summary

Protein kinases are involved in most aspects of cell life due to their role in signal transduction. Dysregulated kinases can cause severe diseases such as cancer, inflammation, and neurodegeneration, which has made them a frequent target in drug discovery for the last decades [@Cohen:2021].
The immense research on kinases has led to an increasing amount of kinase resources [@Kooistra:2017].
Among them is the KLIFS database, which focuses on storing and analyzing structural data on kinases and interacting drugs and other small molecules [@Kanev:2021].
The OpenCADD-KLIFS Python module offers a convenient integration of the KLIFS data into workflows to facilitate computational kinase research.

[OpenCADD-KLIFS](https://opencadd.readthedocs.io/en/latest/databases_klifs.html) (``opencadd.databases.klifs``) is a part of the [OpenCADD](https://opencadd.readthedocs.io/) package, a collection of Python modules for structural cheminformatics.

# Statement of need

The KLIFS resource [@Kanev:2021] contains information about kinases, structures, ligands, interaction fingerprints, and bioactivities. 
KLIFS thereby focuses especially on the ATP binding site, defined as a set of 85 residues and aligned across all structures using a multiple sequence alignment (MSA) [@vanLinden:2014].
Fetching, filtering, and integrating the KLIFS content on a larger scale into Python-based pipelines is currently not straight-forward, especially for users without a background in online queries. Effortless switching between data queries from a _local_ KLIFS download and the _remote_ KLIFS database is not possible.

OpenCADD-KLIFS is aimed at current and future users of the KLIFS database who seek to 
integrate kinase resources into Python-based research projects.
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or remotely from the KLIFS webserver. 
The presented module provides identical APIs for the remote and local queries and streamlines all output into 
standardized Pandas DataFrames [@pandas] to allow for easy and quick downstream data analyses (\autoref{fig:opencadd_klifs_toc}). This Pandas-focused setup is ideal to work with in Jupyter notebooks [@Kluyver:2016]. 

![OpenCADD-KLIFS fetches KLIFS data [@Kanev:2021] offline from a KLIFS download or online from the KLIFS database and formats the output as user-friendly Pandas DataFrames [@pandas].\label{fig:opencadd_klifs_toc}](opencadd_klifs_toc.png)

# State of the field

The KLIFS database is unique in the structure-based kinase field in terms of integrating and annotating different data resources in a kinase- and pocket-focused manner. Kinases, structures, and ligands have unique identifiers in KLIFS, which makes it possible to fetch and filter cross-referenced information for a query kinase, structure, or ligand.

- Kinase structures are fetched from the PDB, split by chains and alternate models, annotated with the KLIFS pocket of 85 residues, and aligned across the fully structurally covered kinome.
- Kinase-ligand interactions seen in experimental structures are annotated for the 85 pocket residues in the form of the KLIFS interaction fingerprint (KLIFS IFP).
- Bioactivity data measured against kinases are fetched from ChEMBL [@Mendez:2018] and linked to kinases, structures, and ligands available in KLIFS.
- Kinase inhibitor metadata are fetched from the PKIDB [@Carles:2018] and linked to co-crystallized ligands available in KLIFS.

The KLIFS data integrations and annotations can be accessed in different ways, which are all open-sourced:

- Manually via the [KLIFS website](https://klifs.net/) interface: This mode is preferable when searching for information on a specific structure or smaller set of structures.
- Automated via the [KLIFS KNIME](https://github.com/3D-e-Chem/knime-klifs) nodes [@McGuire:2017; @Kooistra:2018]: This mode is extremely useful if the users' projects are embedded in KNIME; programming is not needed.
- Programmatically using the REST API and KLIFS OpenAPI specifications: This mode is needed for users who seek to perform larger scale queries or integrate different queries into programmatic workflows. In the following, we will discuss this mode in context of Python-based projects and explain how OpenCADD-KLIFS improves the user experience.

The KLIFS database offers standardized URL schemes (REST API), which allows users to query data by defined URLs, using e.g. the Python package requests [@requests]. Instead of writing customized scripts to generate such KLIFS URLs, the KLIFS OpenAPI specifications &mdash; a document that defines the KLIFS REST API scheme &mdash; can be used to generate a Python client, using e.g. the Python package bravado [@bravado]. This client offers a Python API to send requests and receive responses.
This setup is already extremely useful, however, it has a few drawbacks: the setup is technical, the output is not easily readable for humans and not ready for immediate down-stream integrations &mdash; requiring similar but not identical reformatting functions for different query results &mdash;, and switching from remote requests to local KLIFS download queries is not possible. Facilitating and streamlining these tasks is the purpose of OpenCADD-KLIFS as discussed in more detail in the next section.

# Key Features

The KLIFS database offers a REST API compliant with the OpenAPI specification [@klifsswagger]. Our module OpenCADD-KLIFS uses bravado to dynamically generate a Python client based on the OpenAPI definitions and adds wrappers to enable the following functionalities:

- A session is set up, which allows access to various KLIFS *data sources* by different *identifiers* with the API ``session.data_source.by_identifier``. *Data sources* currently include kinases, structures and annotated conformations, modified residues, pockets, ligands, drugs, and bioactivities; *identifiers* refer to kinase names, PDB IDs, KLIFS IDs, and more. For example, ``session.structures.by_kinase_name`` fetches information on all structures for a query kinase.
- The same API is used for local and remote sessions, i.e. interacting with data from a KLIFS download folder and from the KLIFS website, respectively.
- The returned data follows the same schema regardless of the session type (local/remote); all results obtained with bravado are formatted as Pandas DataFrames with standardized column names, data types, and handling of missing data.
- Files with the structural 3D coordinates deposited on KLIFS include full complexes or selections such as proteins, pockets, ligands, and more. These files can be downloaded to disc or loaded via biopandas [@Raschka:2017] or RDKit [@rdkit].

OpenCADD-KLIFS is especially convenient whenever users are interested in multiple or more complex queries such as "fetching all structures for the kinase EGFR in the DFG-in conformation" or "fetching the measured bioactivity profiles for all ligands that are structurally resolved in complex with EGFR". Formatting the output as DataFrames facilitates subsequent filtering steps and DataFrame merges in case multiple KLIFS datasets need to be combined.

OpenCADD-KLIFS is currently used in several projects from the Volkamer Lab [@volkamerlab] including TeachOpenCADD [@teachopencadd], OpenCADD-pocket [@opencadd_pocket], KiSSim [@kissim], KinoML [@kinoml], and PLIPify [@plipify].
For example, OpenCADD-KLIFS is applied in a [TeachOpenCADD tutorial](https://projects.volkamerlab.org/teachopencadd/talktorials/T012_query_klifs.html) to demonstrate how to fetch all kinase-ligand interaction profiles for all available EGFR kinase structures to visualize the per-residue interaction types and frequencies with only a few lines of code.

# Acknowledgements

We thank the whole KLIFS team for providing such a great kinase resource with an easy-to-use API and especially Albert Kooistra for his help with questions and wishes regarding the KLIFS database. 
We thank David Schaller for his feedback on the OpenCADD-KLIFS module.
We acknowledge the contributors involved in software programs and packages used by OpenCADD-KLIFS, such as bravado, RDKit, Pandas, Jupyter, and Pytest, and Sphinx. 

# References