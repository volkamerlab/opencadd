#!/bin/bash
uname -a
df -h
ulimit -a

. devtools/scripts/initialize_conda.sh
conda info
python devtools/scripts/create_conda_env.py -n=test -p=${{ matrix.python-version }} devtools/conda-envs/test_env.yaml