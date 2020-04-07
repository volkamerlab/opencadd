For developers
==============

Using Atomium
-------------

**Opening a pdb file:**

.. code-block:: python

    pdb1 = atomium.open('../1LOL.pdb')
    pdb2 = atomium.fetch('5XME.pdb')

With "open" you can open a local pdb file.

With "fetch" you can load a pdb file from the RCSB.


**How to select a subset of the protein (all isoleucines, all isoleucines and all valines):**

.. code-block:: python

    model = pdb1.model
    model.residues(name = "ISO")
    model.residues(name__regex = "ISO|VAL")


**Accessing all alpha carbons in the protein:**

.. code-block:: python

    model.atoms(name="CA")


**Removing all waters:**

.. code-block:: python

    model.dehydrate()


**Removing all ligands:**

| There is no function for this in atomium.
| There are to options how to handle this.

**First** we write out one function an bound it to our model-object

.. code-block:: python

    import types

    def remove_ligands(self):
        self._ligands = atomium.base.StructureSet()

    atomium.Model.remove_ligands = types.MethodType(remove_ligands, atomium.Model)

    model.remove_ligands()


**Second**, we can just do it manually

.. code-block:: python

        model._ligands = atomium.base.StructureSet()


**If the pdb file provides a sequence, you can access it with (here for chain A):**

.. code-block:: python

    model.chain("A").sequence

**Computing RMSD (works only, if the structures have the same amount of atoms):**

.. code-block:: python

    model1.rmsd_with(model2)

Alignment API
-------------

This sections defines our input and output parameters and their names.
They should be consistent over all modules.

**How does the API work?**

The first 2 to ``n`` arguments are the structures to be superposed.
If there are more than 2 structures, the first one will be considered the target
structure and every other structure will be superposed on this creating multiple
structural alignments.
If the user does not give any other arguments, the method will be Theseus.
The next argument is the method of choice. If available the user can give some
more method-options.

The structures are automaticaly transformed to `atomium.models
<https://atomium.samireland.com/api/structures.html#atomium.structures.Model>`_.
This means the input for all other alignment functions should be `atomium.models
<https://atomium.samireland.com/api/structures.html#atomium.structures.Model>`_.

**Example:**

.. code-block:: bash

    structuralalignment 1lol 5XME --method mmligner --method-options something


**The aligment function should have the follwing parameters:**

- Input parameters:
    - list of `atomium Models
      <https://atomium.samireland.com/api/structures.html#atomium.structures.Model>`_
    - optional kwargs

- Returns:
    - superposed atomium Models/Molecules
    - scores (RMSD)
    - metadata

How the package will be structured
----------------------------------

Contributing to structuralalignment
-----------------------------------

First of all, thank you for considering making this project better.

You can help to make structuralalignment better by raising an issue to report a bug or suggest a
feature which is not implemented yet.
You can also create a new feature and add it via a Pull Request.

How Pull Requests work
----------------------

1. Fork the ``structuralalignment`` repository into your profile. Now you have your own copy of the repository.
2. Clone the repository. Run ``clone https://github.com/<your_username>/structuralalignment`` in your terminal.
3. Create a new branch to work on: ``git checkout -b meaningful-branch-name``.
4. Make your changes.
5. Add, commit and push your change to your forked repository.
6. Click on ``Create pull request`` Button.
7. Follow the template given there. Be aware of the requirements_.

The ``structuralalignment`` team will review the pull request and if there are no flaws it will be merged
to the ``master`` branch. If there are some flaws, the team will request you to fix them.

.. _requirements:

**************************************
What should each Pull Request contain?
**************************************

* Documentation (should be in ``rst`` format)
* Tests_ with Pytest
* Examples for the use of the new functions
* Little benchmark
* If 3rd party dependencies have been added, rationale behind that decision
* It should be formatted by black with ``-l 99``
* Short summary of the changes that were made

A template is following soon.

Formatting with black
---------------------

**1. Option:**

* Apply "black -l 99" on your code before committing::

        $> black -l 99 structuralalignment/


**2. Option:**

* Configuring your IDE to that effect
* Example in VS Code:

    * go to the settings
    * search for "python formatting"
    * set "Python › Formatting: Provider" to black
    * add "-l 99" to "Python › Formatting: Black Args"
    * activate "Editor: Format On Save"


.. _Tests:

How to add a new test
---------------------

- write a unit test for the new (or changed) function with `pytest
  <https://docs.pytest.org/en/latest/>`_.
- add new dependencies to ``test_env.yaml``


Steps made by the Github Actions Workflow
-----------------------------------------

The actions are executed automatically for every Pull Request submitted,
and for every commit pushed to ``master``. Steps run in Ubuntu and MacOS are:

* Report additional information about the test-build.
* Fixing conda in MacOS (to get the project)
* Creating the environment and getting all necessary dependencies.
* Installing the package in this environment.
* Running the tests.

The formating check is done in ubuntu.

* Checkout the code.
* Installing the linter (pylint) and the formatter (black).
* Running pylint (using  configuration at ``.pylintrc``
* Running ``black -l 99`` in check mode
