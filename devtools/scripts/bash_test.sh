#!/bin/bash
uname -a
df -h
ulimit -a

. initialize_conda.sh
conda activate
conda info
python create_conda_env.py -n=test -p=${{ matrix.python-version }} ../conda-envs/test_env.yaml
conda activate test
pytest -v --cov=superposer --cov-report=xml --color=yes superposer/tests/
conda install -y pylint black
pylint superposer/
black --check -l 99 superposer/