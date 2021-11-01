OpenCADD-KLIFS
================

Once you have installed the ``opencadd`` package, you will have access (among others) 
to the ``opencadd.databases.klifs`` module (OpenCADD-KLIFS).
In case you wish to install only the dependencies relevant to OpenCADD-KLIFS, please follow the installation instructions `here <https://opencadd.readthedocs.io/en/latest/installing_opencadd_klifs.html>`_.

OpenCADD-KLIFS offers a simple API to interact with data from KLIFS remotely and locally.

Find a detailed tutorial at the `TeachOpenCADD platform <https://projects.volkamerlab.org/teachopencadd/talktorials/T012_query_klifs.html>`_ on the KLIFS database and on how to apply the module OpenCADD-KLIFS to an example research question.

What is KLIFS and who created it?
---------------------------------

  KLIFS is a kinase database that dissects experimental structures of catalytic kinase domains and the way kinase inhibitors interact with them. The KLIFS structural alignment enables the comparison of all structures and ligands to each other. Moreover, the KLIFS residue numbering scheme capturing the catalytic cleft with 85 residues enables the comparison of the interaction patterns of kinase-inhibitors, for example, to identify crucial interactions determining kinase-inhibitor selectivity.

- KLIFS database: https://klifs.net (official), https://dev.klifs.net/ (developmental)
- KLIFS online service: https://klifs.net/swagger (official), https://dev.klifs.net/swagger_v2 (developmental, used here) 
- KLIFS citation: `Nucleic Acids Res. (2021), 49, D1, D562–D569 <https://academic.oup.com/nar/article/49/D1/D562/5934416>`_

What does ``opencadd.databases.klifs`` offer?
---------------------------------------------

This module allows you to access KLIFS data such as information about 
kinases, structures, structural conformations, modified residues, ligands, drugs, interaction fingerprints, and bioactivities.

On the one hand, you can query the KLIFS webserver directly. 
On the other hand, you can query your local KLIFS download.
We provide identical APIs for the remote and local queries and streamline all output into standardized ``pandas`` DataFrames for easy and quick downstream data analyses.

Work with KLIFS data from KLIFS server (remotely)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``opencadd.databases.klifs.remote`` submodule offers you to access KLIFS data from the KLIFS server.

Our API relies on the REST API and OpenAPI (formerly Swagger API) specification at https://dev.klifs.net/swagger_v2/ to dynamically generate a Python client with ``bravado``.

Example for ``opencadd``'s API to access remote data:

.. code-block:: python

    from opencadd.databases.klifs import setup_remote

    # Set up remote session
    session = setup_remote()

    # Get all kinases that are available remotely
    session.kinases.all_kinases()

    # Get kinases by kinase name
    session.kinases.by_kinase_name(["EGFR", "BRAF"])

Work with KLIFS data from disc (locally)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``opencadd.databases.klifs.local`` submodule offers you to access KLIFS data from your KLIFS download. 
In order to make use of the module's functionality, you need a KLIFS download folder ``KLIFS_download`` with the following structure (files downloaded from `KLIFS <from https://klifs.net>`_):

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

Example for ``opencadd``'s API to access local data:

.. code-block:: python

    from opencadd.databases.klifs import setup_local

    # Set up local session
    session = setup_local("../../opencadd/tests/databases/data/KLIFS_download")

    # Get all kinases that are available locally
    session.kinases.all_kinases()

    # Get kinases by kinase name
    session.kinases.by_kinase_name(["EGFR", "BRAF"])


How is ``opencadd.databases.klifs`` structured?
----------------------------------------------------------

The module's structure looks like this, using the same API for both modules ``local`` and ``remote`` whenever possible:

.. code-block:: console 

    opencadd/ 
        └── databases/
            └── klifs/
                ├── api.py         # Defines the main API for local and remote sessions.
                ├── session.py     # Defines a KLIFS session.
                ├── core.py        # Defines the parent classes used in the local and remote modules.
                ├── local.py       # Defines the API for local queries.
                ├── remote.py      # Defines the API for remote queries.
                ├── schema.py      # Defines the schema for class method return values.
                ├── fields.py      # Defines the different KLIFS data fields and their names/dtypes in ``opencadd``.
                ├── utils.py       # Defines utility functions.
                └── exceptions.py  # Defines exceptions.

This structure mirrors the KLIFS OpenAPI structure in the following way to access different kinds of information both remotely and locally:

- ``kinases``  

  - Get information about kinases (groups, families, names).  
  - In KLIFS OpenAPI called ``Information``: https://dev.klifs.net/swagger_v2/#/Information

- ``ligands``  

  - Get ligand information.  
  - In KLIFS OpenAPI called ``Ligands``: https://dev.klifs.net/swagger_v2/#/Ligands

- ``structures``

  - Get structure information.  
  - In KLIFS OpenAPI called ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures  

- ``bioactivities``  

  - Get bioactivity information.  
  - In KLIFS OpenAPI part of ``Ligands``: https://dev.klifs.net/swagger_v2/#/Ligands  

- ``interactions``  

  - Get interaction information.  
  - In KLIFS OpenAPI called ``Interactions``: https://dev.klifs.net/swagger_v2/#/Interactions  

- ``pocket``  

  - Get interaction information.  
  - In KLIFS OpenAPI part of ``Interactions``: https://dev.klifs.net/swagger_v2/#/Interactions/get_interactions_match_residues 

- ``coordinates``  

  - Get structural data (structure coordinates).
  - In KLIFS OpenAPI part of ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures 

- ``conformations``

  - Get information on structure conformations.
  - In KLIFS OpenAPI part of ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures/get_structure_conformation

- ``modified_residues``

  - Get information on residue modifications in structures.
  - In KLIFS OpenAPI part of ``Structures``: https://dev.klifs.net/swagger_v2/#/Structures/get_structure_modified_residues


