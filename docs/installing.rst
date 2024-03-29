Installing
==========

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 

    If you installed ``mamba`` into an existing ``conda`` installation, also make sure that the ``conda-forge`` channel is configured by running ``conda config --add channels conda-forge``.


Install from the conda package
------------------------------

1. Create a new conda environment called ``opencadd`` with the ``opencadd`` package and all its dependencies installed::

    mamba create -n opencadd opencadd

2. Activate the new conda environment::

    conda activate opencadd

.. 3. Test that your installation works::

    superposer -h


Install from the latest development snapshot
--------------------------------------------

Install the latest development snapshot from the `GitHub repository's master branch <https://github.com/volkamerlab/opencadd>`_.


1. Create a new conda environment called ``opencadd``::

    mamba env create -f https://raw.githubusercontent.com/volkamerlab/opencadd/master/devtools/conda-envs/user_env.yaml

2. Activate the new conda environment::

    conda activate opencadd

3. Install ``opencadd`` package via pip::

    pip install https://github.com/volkamerlab/opencadd/archive/master.tar.gz

.. 4. Test that your installation works::

    superposer -h


Development version
-------------------

To install the development version of OpenCADD, you can run::

    # Install development environment (incl. packages for testing and documentation)
    mamba env create -f https://raw.githubusercontent.com/volkamerlab/opencadd/master/devtools/conda-envs/test_env.yaml -n opencadd-dev
    conda activate opencadd-dev
    
    # Download repository and install `opencadd` package in editable mode
    git clone git@github.com:volkamerlab/opencadd.git
    pip install -e opencadd

    # Let's change into the repo's folder
    cd opencadd
    
    # Run tests like this
    pytest -v opencadd/tests/

    # Before you push changes to GitHub, lint and style-check your code
    pylint opencadd
    black --check -l 99 opencadd

    # Check how the documentation renders
    cd docs
    make html

Note: If you add new dependencies to ``opencadd``, please update the 
`test <https://github.com/volkamerlab/opencadd/blob/master/devtools/conda-envs/test_env.yaml>`_ and 
`user <https://github.com/volkamerlab/opencadd/blob/master/devtools/conda-envs/user_env.yaml>`_ environment, 
and leave a note in the 
`dependencies <https://github.com/volkamerlab/opencadd/blob/master/docs/installing.rst#dependencies>`_ section.


Dependencies
------------

``opencadd`` is supported for at least Python 3.7 and needs the following modules: 

+---------------------+---------------------------+--------------------+--+--+
| Package             | Minimal supported version | Used by submodules |  |  |
+=====================+===========================+====================+==+==+
| ``pandas``          | 1.3.4                     | 1, 2, 4            |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``biopandas``       | 0.2.9                     | 1, 2               |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``bravado``         | 11.0.3                    | 1                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``rdkit``           | 2021.09.2                 | 1                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``tqdm``            | 4.62.3                    | 1                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``biopython``       | = 1.77                    | 2                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``mdanalysis``      | 2.0.0                     | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``biotite``         | 0.31.0                    | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``mmligner``        | 1.0.2                     | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``muscle``          | 3.8.1551                  | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``theseus``         | 3.3.0                     | 3                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``matplotlib-base`` | 3.5.0                     | 4                  |  |  |
+---------------------+---------------------------+--------------------+--+--+
| ``nglview``         | 3.0.3                     | 4                  |  |  |
+---------------------+---------------------------+--------------------+--+--+


Some packages are only needed for a subset of the following modules: [1] ``opencadd.databases.klifs``, 
[2] ``opencadd.io``, 
[3] ``opencadd.structure.superposition``, 
[4] ``opencadd.structure.pocket``

This list of minimal supported versions is based on `this CI run <https://github.com/volkamerlab/opencadd/runs/4462667598?check_suite_focus=true#step:6:42>`_.