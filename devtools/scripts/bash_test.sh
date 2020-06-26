#!/bin/bash
uname -a
df -h
ulimit -a

. initialize_conda.sh
conda activate
conda info
python create_conda_env.py -n=test -p=${{ matrix.python-version }} ../conda-envs/test_env.yaml
conda activate test
pytest -v --cov=opencadd --cov-report=xml --color=yes opencadd/tests/
conda install -y pylint black
pylint opencadd/
black --check -l 99 opencadd/