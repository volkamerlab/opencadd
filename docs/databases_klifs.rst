Databases: KLIFS
================

Once you have installed the package, you will have access (among others) 
to the ``opencadd.databases.klifs`` module.

This module offers a simple API to interact with data from KLIFS remotely and locally.


What is KLIFS and who created it?
---------------------------------

"Kinase-Ligand Interaction Fingerprints and Structures database (KLIFS), developed at the Division of Medicinal Chemistry - VU University Amsterdam, is a database that revolves around the protein structure of catalytic kinase domains and the way kinase inhibitors can interact with them."

- KLIFS database: https://klifs.vu-compmedchem.nl 
- KLIFS online service: https://klifs.vu-compmedchem.nl/swagger 
- KLIFS citation 

  - Description of the database itself: `J. Med. Chem. (2014), 57, 2, 249-277 <https://pubs.acs.org/doi/abs/10.1021/jm400378w>`_ 
  - Description of the online service: `Nucleic Acids Res. (2016), 44, 6, D365–D371 <https://academic.oup.com/nar/article/44/D1/D365/2502606>`_ 


What does ``opencadd.databases.klifs`` offer?
---------------------------------------------

Work with KLIFS data from KLIFS server (remotely)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``opencadd.databases.klifs.remote`` submodule offers you to access KLIFS data from the KLIFS server.

This module uses the official KLIFS API: https://klifs.vu-compmedchem.nl/swagger.

.. code-block:: python

    from opencadd.databases.klifs import setup_remote

    # Set up remote session
    remote = setup_remote()

    # Get all kinases that are available remotely
    remote.kinases.all_kinases()

    # Get kinases by kinase name
    remote.kinases.from_kinase_names(["EGFR", "BRAF"])

Work with KLIFS data from disc (locally)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``opencadd.databases.klifs.local`` submodule offers you to access KLIFS data from the KLIFS server. In order to make use of the module's functionality, you need a KLIFS download folder ``KLIFS_download`` with the following structure (files downloaded from `KLIFS <from https://klifs.vu-compmedchem.nl>`_):

.. code-block:: console 

    └── KLIFS_download 
        ├── KLIFS_export.csv           # Metadata file 
        ├── overview.csv               # Metadata file 
        └── HUMAN     	               # Species name 
            ├── AAK1                   # Kinase name 
            │   ├── 4wsq_altA_chainA   # PDB ID, alternate model ID, chain ID 
            │   │   ├── complex.mol2 
            │   │   ├── ligand.mol2 
            │   │   ├── protein.mol2 
            │   │   ├── pocket.mol2 
            │   │   └── water.mol2 
            │   └── ... 
            └── ... 

.. code-block:: python

    from opencadd.databases.klifs import setup_local

    # Set up local session
    local = setup_local("../../opencadd/tests/databases/data/KLIFS_download")

    # Get all kinases that are available locally
    local.kinases.all_kinases()

    # Get kinases by kinase name
    local.kinases.from_kinase_names(["EGFR", "BRAF"])


How is ``opencadd.databases.klifs`` structured?
----------------------------------------------------------

The module's structure looks like this, trying to use the same API for both modules ``local`` and ``remote`` whenever possible:

.. code-block:: console 

    opencadd/ 
        └── databases/
            └── klifs/
                ├── api.py
                ├── core.py
                ├── local.py
                ├── remote.py
                ├── schema.py
                └── utils.py

This structure mirrors the KLIFS Swagger API structure in the following way to access different kinds of information both remotely and locally:

- ``kinases``  

  - Get information about kinases (groups, families, names).  
  - In KLIFS swagger API called ``Information``.  

- ``ligands``  

  - Get ligand information.  
  - In KLIFS swagger API called ``Ligands``.  

- ``structures``

  - Get structure information.  
  - In KLIFS swagger API called ``Structures``.  

- ``bioactivities``  

  - Get bioactivity information.  
  - In KLIFS swagger API part of ``Ligands``.  

- ``interactions``  

  - Get interaction information.  
  - In KLIFS swagger API called ``Interactions``.  

- ``pocket``  

  - Get interaction information.  
  - In KLIFS swagger API part of ``Interactions``.  

- ``coordinates``  

  - Get structural data (structure coordinates).
  - In KLIFS swagger API part of ``Structures``.  


