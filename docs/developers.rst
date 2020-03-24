For the developers
==================

Using Atomium
-------------


Opening a pdb file::

    pdb1 = atomium.open('../1LOL.pdb')
    pdb2 = atomium.fetch('5XME.pdb')
With "open" you can open a local pdb file.

With "fetch" you can load a pdb file from the RCSB.


How to select a subset of the protein (all isoleucines, all isoleucines and all valines)::

    model = pdb1.model
    model.residues(name = "ISO")
    model.residues(name__regex = "ISO|VAL")


Accessing all alpha carbons in the protein::

    model.atoms(name="CA")


Removing all waters::

    model.dehydrate()


Removing all ligands::

There is no function for this in atomium.
There are to options how to handle this.
First we write out one function an bound it to our model-object::
    import types

    def remove_ligands(self):
        self._ligands = atomium.base.StructureSet()

    atomium.Model.remove_ligands = types.MethodType(remove_ligands, atomium.Model)

    model.remove_ligands()


Second, we can just do it manually::

        model._ligands = atomium.base.StructureSet()


If the pdb file provides a sequence, you can access it with (here for chain A)::

    model.chain("A").sequence

Computing RMSD (works only, if the structures have the same amount of atoms)::

    model1.rmsd_with(model2)

Alignment API
-------------

This sections defines our input and output parameters and their names.
They should be consistent over all modules.

The aligment function should have the follwing parameters:

- Input parameters:
    - atomium Models [https://atomium.samireland.com/api/structures.html]
    - optional kwargs

- Returns:
    - superposed structures using atomium Models/Molecules
    - scores (RMSD)

How the package will be structured
----------------------------------

What should each Pull Request contain
---------------------------------------

* Documantation
* Tests with Pytest
* Examples for the use of the new functions
* Little benchmark
* Explaination of 3rd party dependencies

A template is following soon.

How to add a new test
---------------------

Coming soon.