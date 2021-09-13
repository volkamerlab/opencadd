Installing
==========

.. note::

    We are assuming you have a working ``mamba`` installation in your computer. 
    If this is not the case, please refer to their `official documentation <https://mamba.readthedocs.io/en/latest/installation.html#mamba>`_. 


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

    mamba env create -f https://raw.githubusercontent.com/volkamerlab/opencadd/master/devtools/conda-envs/test_env.yaml -n opencadd

2. Activate the new conda environment::

    conda activate opencadd

3. Install ``opencadd`` package via pip::

    pip install https://github.com/volkamerlab/opencadd/archive/master.tar.gz

.. 4. Test that your installation works::

    superposer -h
