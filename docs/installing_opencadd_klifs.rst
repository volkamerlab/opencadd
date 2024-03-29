Installing OpenCADD-KLIFS only
==============================

In case you would like to install the dependencies for the OpenCADD-KLIFS module only, please follow these instructions.

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 

    If you installed ``mamba`` into an existing ``conda`` installation, also make sure that the ``conda-forge`` channel is configured by running ``conda config --add channels conda-forge``.


1. Create a new conda environment called ``opencadd-klifs`` only with the dependencies needed for the OpenCADD-KLIFS module::

    mamba create -n opencadd-klifs bravado pandas tqdm rdkit biopandas

2. Activate the new conda environment::

    conda activate opencadd-klifs

3. Install the ``opencadd`` package without any dependencies (all OpenCADD-KLIFS-relevant dependencies have been installed in step 1 already)::

    mamba install opencadd --no-deps

   If you are planning on working with Jupyter notebooks, install JupyterLab and IPyWidgets::

    mamba install jupyterlab ipywidgets
