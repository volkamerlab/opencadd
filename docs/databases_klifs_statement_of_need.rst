Statement of need
================= 

The KLIFS resource [Kanev_2021]_ contains information about kinases, structures, ligands, 
interaction fingerprints, and bioactivities. 
KLIFS thereby focuses especially on the ATP binding site, defined as a set of 85 residues and 
aligned across all structures using a multiple sequence alignment [vanLinden_2014]_.
Fetching, filtering, and integrating the KLIFS content on a larger scale into Python-based 
pipelines is currently not straight-forward, especially for users without a background in 
online queries. 
Furthermore, switching between data queries from a *local* KLIFS download and 
the *remote* KLIFS database is not readily possible.

OpenCADD-KLIFS is aimed at current and future users of the KLIFS database who seek to 
integrate kinase resources into Python-based research projects.
With OpenCADD-KLIFS, KLIFS data can be queried either locally from a KLIFS download or 
remotely from the KLIFS webserver. 
The presented module provides identical APIs for the remote and local queries and 
streamlines all output into standardized Pandas DataFrames 
`Pandas <https://doi.org/10.5281/zenodo.5574486>`_  to allow for easy and quick 
downstream data analyses (Figure 1). 
This Pandas-focused setup is ideal if you work with Jupyter notebooks [Kluyver_2016]_.

.. raw:: html

   <p align="center">
   <img src="_static/opencadd_klifs_toc.png" alt="OpenCADD-KLIFS" width="600"/>
   </p>

*Figure 1*: OpenCADD-KLIFS fetches KLIFS data offline from a KLIFS download or 
online from the KLIFS database and formats the output as user-friendly Pandas DataFrames.

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