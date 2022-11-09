import os
from typing import NoReturn
from pathlib import Path


def pdb_to_pdbqt(
        filepath_input_pdb: Path,
        filepath_output_pdbqt: Path,
        filepath_mgl_python: Path = "/home/michele/Software/MGLTools-1.5.6/bin/pythonsh",
        filepath_mgl_autodock: Path = "/home/michele/Software/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py",
) -> NoReturn:
    """
    Convert a PDB file to a PDBQT file, using AutoDock script in MGLTools.
    The PDBQT file format is the structure format designed by AutoDock, and is the required input
    file format in most AutoDock functionalities.

    Parameters
    ----------
    filepath_input_pdb : pathlib.Path
        Path to the target PDB file.
    filepath_output_pdbqt : pathlib.Path
        Path for storing the output PDBQT file.
    filepath_mgl_python : pathlib.Path
        Path to the Python executable file of the MGLTools program.
    filepath_mgl_autodock : pathlib.Path
        Path to the `prepare_receptor4.py` module of the AutoDock Tools in MGLTools program.

    Returns
    -------
    None
        The converted PDBQT file will be generated in the given location.
    """
    #TODO: This function uses an outdated version of MGLTools, which requires separate installation
    # and manual input of MGLTools paths. Update to a more convenient option!
    os.system(
        f"{filepath_mgl_python} {filepath_mgl_autodock} "
        f"-r {filepath_input_pdb} -o {filepath_output_pdbqt.with_suffix}"
    )
    return
