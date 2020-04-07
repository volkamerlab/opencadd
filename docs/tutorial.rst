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



best pactices
-------------

coming soon.
