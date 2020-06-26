Getting Started
===============

Once you have installed the package, you will have access to a ``opencadd`` executable. Some quick examples:

Two structures superposed with mmligner without any additional options::

    opencadd 1lol 5XME
    opencadd 1lol 5XME --method theseus

Input structures can be PDB identifiers, or paths. to local files. If more than two structures are provided,
the first will be considered the target structure (it will remain static). The other structure will be superposed
on top of this one, creating multiple structure alignments.

``opencadd`` will create several PDB files on disk for you to check with your favourite visualization tool.
