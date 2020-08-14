Databases: KLIFS
================

Once you have installed the package, you will have access to the ``opencadd.databases.klifs`` module (among other utilities!), which offers functions to interact remotely and locally with data from KLIFS.


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

    import opencadd.databases.klifs as klifs

    # Get all kinases that are available remotely
    klifs.remote.kinases.kinase_names()

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

**Note**: In order to use the ``opencadd.databases.klifs.local`` module, it is necessary to initialize a metadata file (**the API needs to be simplified here...**).

.. code-block:: python

    import opencadd.databases.klifs as klifs

    # Initialize KLIFS metadata
    klifs_overview_path = Path("/path/to/overview.csv")
    klifs_export_path = Path("/path/to/KLIFS_export.csv")
    klifs.local.initialize.from_files(klifs_overview_path, klifs_export_path)

    # Load KLIFS metadata
    klifs_metadata_path = Path("path/to/klifs_metadata.csv")
    klifs_metadata = pd.read_csv(klifs_metadata_path)

    # Get all kinases that are available locally
    klifs.remote.kinases.kinase_names(klifs_metadata)


How is ``opencadd.databases.klifs`` structured?
----------------------------------------------------------

The module's structure looks like this, trying to use the same API for both modules ``local`` and ``remote`` whenever possible:

.. code-block:: console 

    opencadd/ 
        └── databases/
            └── klifs/
                ├── local/
                │   ├── coordinates.py
                │   ├── initialize.py
                │   ├── interactions.py
                │   ├── kinases.py
                │   ├── ligands.py
                │   └── structures.py
                ├── remote/
                │   ├── coordinates.py
                │   ├── interactions.py
                │   ├── kinases.py
                │   ├── ligands.py
                │   └── structures.py
                ├── klifs_client.py
                └── utils.py

This structure mirrors the KLIFS Swagger API structure in the following way to access different kinds of information both remotely and locally:

- ``kinases`` (in KLIFS called ``information``): Get information about kinases (groups, families, names)
- ``interactions``: Get interaction fingerprint via structure_ID
- ``ligands``: Get ligand information
- ``structures``: Get structure information
- ``coordinates`` (in KLIFS part of ``structures``): Get structural data (structure coordinates)

