Structure: Pocket
=================

Once you have installed the package, you will have access (among others) 
to the ``opencadd.structure.pocket`` module.

This module offers a simple API to define and visualize 
(``nglview``) subpockets and regions within a pocket.

``nglview``: http://nglviewer.org/nglview/latest/

Pocket
------

The ``Pocket`` class is the main API to work with. A ``Pocket`` object is initialized from a 
protein structure file, a name, pocket residue PDB IDs and optionally pocket residue labels.
For instance, pocket residue labels can be used to label residues with alignment IDs, 
such as the KLIFS residue ID (pocket residues 1-85 are aligned across all kinases).


.. code-block:: python

    from opencadd.structure.pocket import Pocket

    pocket = Pocket.from_file(
        filepath="protein.mol2",  # Dummy file path
        pocket_residue_pdb_ids=[50, 51, ..., 198], 
        name="example kinase", 
        pocket_residue_labels=[1, 2, ..., 85]
    )


Subpockets
----------

Currently, the main focus of this module is the definition of subpockets within a pocket.
Once a ``Pocket`` object is initialized, subpockets can be added one-by-one.
The user can define each subpocket by so-called anchor residues. The centroid of all anchor
residues' CA atoms is the subpocket center and will be visualized as a subpocket sphere 
(together with the anchor residues' CA atoms as small spheres). 
Anchor residues need to be selected manually, in a way that their centroid will cover the subpocket 
center as intended.

As an example, we add here a kinase subpocket called ``"hinge"``, which will be 
visualized in magenta and whose center is calculated based on the CA atoms of residues 73, 128, and
193 (residue PDB IDs). These residue PDB IDs correspond to the KLIFS alignment IDs 16, 47, and 80.


.. code-block:: python

    pocket.add_subpocket(
        name="hinge", 
        anchor_residue_pdb_ids=[73, 128, 193], 
        color="magenta", 
        anchor_residue_labels=[16, 47, 80]  # Optionally
    )


Regions
------- 

Usually, it is helpful to highlight the pocket region, which inspired the subpocket choice, 
hence the module allows definitions of important pocket regions. 

In our example, the subpocket ``"hinge"`` is intended to target one of the most important
regions in kinases, the hinge region, where kinase inhibitors form hydrogen bonds with the kinase.
We use the hinge region residues to add a region to our ``Pocket`` object.

.. code-block:: python

    pocket.add_region(
        name="hinge region", 
        residue_pdb_ids=[127, 128, 129], 
        color="magenta", 
        residue_labels=[46, 47, 48]  # Optionally
    )


Visualize the pocket
--------------------

Now we can visualize the pocket using:

.. code-block:: python

    pocket.visualize()


Check out our tutorial to find out more!
