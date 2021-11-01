Installing OpenCADD-KLIFS only
==============================

In case you would like to install the dependencies for the OpenCADD-KLIFS module only, please follow these instructions.

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 


Install from the conda package
------------------------------

1. Create a new conda environment called ``opencadd-klifs`` with the ``opencadd`` package and all its dependencies installed::

    mamba create -n opencadd-klifs bravado pandas tqdm rdkit biopandas

2. Activate the new conda environment::

    conda activate opencadd-klifs

3. Install ``opencadd`` without any dependencies (all ``opencadd-klifs`` relevant dependencies have been installed in step 1)::

    mamba install opencadd --no-deps

   If you are planning on working with Jupyter notebooks, install JupyterLab::

    mamba install jupyterlab ipywidgets
