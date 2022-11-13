import os
from typing import NoReturn, Tuple, Optional, Union
from pathlib import Path

from openbabel import pybel


def pdb_to_pdbqt_autodocktools(
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


def pdb_to_pdbqt_openbabel(
        pdb_filepath: Path,
        pdbqt_filepath: Optional[Path] = None,
        add_hydrogens: bool = True,
        protonate_for_pH: Union[float, None] = 7.4,
        calculate_partial_charges: bool = True,
):
    """
    Convert a PDB file to a PDBQT file, and save it in the given filepath.

    Parameters
    ----------
    pdb_filepath: str or pathlib.Path
        Path to input PDB file.
    pdbqt_filepath: str or pathlib.Path
        Path to output PDBQT file.
    add_hydrogens : bool, Optional, default: True
        Whether to add hydrogen atoms to the structure.
    protonate_for_pH : float | None, Optional, default: 7.4
        pH value to optimize protonation state of the structure. Disabled if `None`.
    calculate_partial_charges : bool, Optional, default: True
        Whether to calculate partial charges for each atom.

    Returns
    -------
    openbabel.pybel.Molecule
        Molecule object of PDB file, modified according to the input.
        The PDBQT file will be stored in the provided path.
    """
    # pybel.readfile() provides an iterator over the Molecules in a file.
    # To access the first (and possibly only) molecule in a file, we use next()
    molecule = next(pybel.readfile("pdb", str(Path(pdb_filepath).with_suffix(".pdb"))))
    if protonate_for_pH:
        molecule.OBMol.CorrectForPH(protonate_for_pH)
    if add_hydrogens:
        molecule.addh()
    if calculate_partial_charges:
        for atom in molecule.atoms:
            atom.OBAtom.GetPartialCharge()
    molecule.write("pdbqt", str(pdbqt_filepath.with_suffix(".pdbqt")), overwrite=True)
    return molecule
