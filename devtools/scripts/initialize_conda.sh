
case ${CI_OS} in
    windows*)
        eval "$(${CONDA}/condabin/conda.bat shell.bash hook)";;
    macOS*)
        eval "$(${CONDA}/condabin/conda shell.bash hook)";;
    *)
        eval "$(conda shell.bash hook)";;
esac

if [ -d ${CONDA}/envs/test ] || [ -d $HOME/.conda/envs/test ]; then
    conda activate test
else
    conda activate
fi