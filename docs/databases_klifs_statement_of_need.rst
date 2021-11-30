Statement of need
=================



`OpenCADD-KLIFS <https://opencadd.readthedocs.io/en/latest/databases_klifs.html>`_ 
(``opencadd.databases.klifs``) is a part of the `OpenCADD <https://opencadd.readthedocs.io/>`_ 
package, a collection of Python modules for structural cheminformatics.
This module offers access to KLIFS data [Kanev_2021]_ such as information about kinases, 
structures, ligands, 
interaction fingerprints, and bioactivities. 
KLIFS thereby focuses especially on the ATP binding site, defined as a set of 85 residues and 
aligned across all structures using a multiple sequence alignment (MSA) [vanLinden_2014]_.
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or remotely 
from the KLIFS webserver. 
The presented module provides identical APIs for the remote and local queries for KLIFS data and 
streamlines all output into 
standardized `Pandas <https://doi.org/10.5281/zenodo.5574486>`_ DataFrames to allow for easy and quick downstream data analyses 
(Figure 1). This Pandas-focused setup is ideal to work with in Jupyter 
notebooks [Kluyver_2016]_. 


.. raw:: html

   <p align="center">
   <img src="_static/opencadd_klifs_toc.png" alt="OpenCADD-KLIFS" width="600"/>
   </p>

*Figure 1*: OpenCADD-KLIFS fetches KLIFS data offline from a KLIFS download or 
online from the KLIFS database and formats the output as user-friendly Pandas DataFrames.

The KLIFS database offers a REST API compliant with the OpenAPI specification 
(`KLIFS OpenAPI <https://dev.klifs.net/swagger_v2/>`_). 
Our module OpenCADD-KLIFS uses `bravado <https://github.com/Yelp/bravado>`_ to dynamically 
generate a Python client based on the OpenAPI definitions and adds wrappers to enable the 
following functionalities:

- A session is set up, which allows access to various KLIFS *data sources* by different 
  *identifiers* with the API ``session.data_source.by_identifier``. *Data sources* currently 
  include kinases, structures and annotated conformations, modified residues, pockets, ligands, 
  drugs, and bioactivities; *identifiers* refer to kinase names, PDB IDs, KLIFS IDs, and more. 
  For example, ``session.structures.by_kinase_name`` fetches information on all structures for a 
  query kinase.
- The same API is used for local and remote sessions.
- The returned data follows the same schema regardless of the session type (local/remote); all 
  results obtained with bravado are formatted as Pandas DataFrames with standardized column names, 
  data types, and handling of missing data.
- Files with the structural 3D coordinates deposited on KLIFS include full complexes or selections 
  such as proteins, pockets, ligands, and more. These files can be downloaded to disc or loaded 
  via biopandas [Raschka_2017]_ or `RDKit <http://www.rdkit.org>`_. 

OpenCADD-KLIFS is especially convenient whenever users are interested in multiple or more 
complex queries such as "fetching all structures for the kinase EGFR in the DFG-in conformation" 
or "fetching the measured bioactivity profiles for all ligands that are structurally resolved in 
complex with EGFR". Formatting the output as DataFrames facilitates subsequent filtering steps 
and DataFrame merges in case multiple KLIFS datasets need to be combined.
OpenCADD-KLIFS is currently used in several projects 
from the `Volkamer Lab <https://volkamerlab.org/>`_ 
including 
`TeachOpenCADD <https://github.com/volkamerlab/teachopencadd>`_, 
`OpenCADD-pocket <https://github.com/volkamerlab/opencadd>`_, 
`KiSSim <https://github.com/volkamerlab/kissim>`_, 
`KinoML <https://github.com/openkinome/kinoml>`_, and 
`PLIPify <https://github.com/volkamerlab/plipify>`_.
For example, OpenCADD-KLIFS is applied in a 
`TeachOpenCADD tutorial <https://projects.volkamerlab.org/teachopencadd/talktorials/T012_query_klifs.html>`_ 
to demonstrate how to fetch all kinase-ligand interaction profiles for all available EGFR kinase 
structures to visualize the per-residue interaction types and frequencies with only a few 
lines of code.

.. [Kanev_2021] Kanev et al., (2021),
   KLIFS: an overhaul after the first 5 years of supporting kinase research,
   Nucleic Acids Research, 
   49(D1), D562–D569, doi:10.1093/nar/gkaa895.
.. [vanLinden_2014] van Linden et al., (2014)
   KLIFS: A Knowledge-Based Structural Database To Navigate Kinase–Ligand 
   Interaction Space, 
   Journal of Medicinal Chemistry, 
   57(2), 249-277, doi:10.1021/jm400378w.
.. [Kluyver_2016] Kluyver et al., (2016),
   Jupyter Notebooks – a publishing format for reproducible computational workflows,
   In Positioning and Power in Academic Publishing: Players, Agents and Agendas. IOS Press. pp. 87-90,
   doi:10.3233/978-1-61499-649-1-87.
.. [Raschka_2017] Raschka, (2017), 
   BioPandas: Working with molecular structures in pandas DataFrames, Journal of Open Source Software, 
   2(14), 279, doi:10.21105/joss.00279.