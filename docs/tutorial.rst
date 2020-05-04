Tutorials
=========

This page shows the user how to use this library.

Important notes
---------------

If the imput is more than 2 structures, the first will be considered the target structure.
Meaning every other structure will be superposed with this one, creating multiple structure alignments.
The sequence length has to be the same.

.. The superposed structures will be saved in a pdb file.

How to use the library
----------------------

After the command :code:`structuralalignment` the user provides the structurenames. The library will
search locally for the file, if a path is given. Otherwise it will search in the `Protein Data Bank
<https://www.rcsb.org/>`_.

Next the user will choose the wanted method with :code:`--method [name of the method]`.

Optionally the user can also configure the method options with :code:`--method-options
[name of the option: new value; name of the option: new value; ...]`.

Simple examples
---------------
Two structures superposed with mmligner without any additional options::

    structuralalignment 1lol 5XME --method mmligner


Two structures superposed with theseus with additional arguments::

    structuralaligment 1lol 5XME --method theseus --method-options something

Two structures superposed with theseus with additional options::

    structuralalignment 1lol 5XME --method matchmaker --method-strategy --method-matrix --method-gap --method-strict --method-selection --method-weights --method-mass_tolerance



best pactices
-------------

MatchMaker:
    Introduction:
        MatchMaker superimposes protein by first creating pairwise sequence alignments, then fitting the aligned residue pairs.
        The standard Needleman-Wunsch and Smith-Waterman algorithms are available for producing global and local sequence alignments
        and MatchMaker can identify the best-matching chains based on alignment scores.
        Alignment scores can include residue similarity, secondary structure information, and gap penalties.
        MatchMaker performs a spatial superposition by minimizing the RMSD.
    Preparation:
        Generating pairwise sequence alignments and matching, i.e., superimposing the structures according to those pairwise alignments
    Alignment:
        Spatially align the group of atoms `mobile` to `reference` by doing a RMSD fit on `select` atoms.
        - A rotation matrix is computed that minimizes the RMSD between
        the coordinates of `mobile.select_atoms(sel1)` and
        `reference.select_atoms(sel2)`; before the rotation, `mobile` is
        translated so that its center of geometry (or center of mass)
        coincides with the one of `reference`.
        - All atoms in :class:`~MDAnalysis.core.universe.Universe` that
        contain `mobile` are shifted and rotated.
    Analysis:
        RMSD before and after spatial alignment


