Structural superposition
========================

.. todo:
    Consider using ``opencadd superpose`` as a subcommand.

Once you have installed the package, you will have access to a ``superposer`` executable (among other utilities!). Some quick examples:

Two structures superposed with mmligner without any additional options::

    superposer 1lol 5XME
    superposer 1lol 5XME --method theseus

Input structures can be PDB identifiers, or paths. to local files. If more than two structures are provided,
the first will be considered the target structure (it will remain static). The other structure will be superposed
on top of this one, creating multiple structure alignments.

``superposer`` will create several PDB files on disk for you to check with your favourite visualization tool.
