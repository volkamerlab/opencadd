Installing
==========

Eventually, we will have a ``conda`` package, but for now you need to create a new environment manually.

1. Install Miniconda for your OS if you don't have it already.
2. Download a `copy of this repository <https://github.com/volkamerlab/opencadd/archive/master.zip>`_.
3. Create new conda environment::

    conda env create -n opencadd -f devtools/conda-envs/test_env.yaml

4. Activate the new environment::

    conda activate opencadd

5. Install the package with ``pip``::

    python -m pip install https://github.com/volkamerlab/opencadd/archive/master.tar.gz

6. Run ``superposer -h`` to test it works.

7. Workaround for MMLigner::
    conda config --add channels conda-forge 
    conda activate base
    conda install conda-build
    conda build devtools/conda-recipes/mmligner/
    (conda activate opencadd)
    conda install -c local mmligner pip
    python -m pip install -e . --no-deps