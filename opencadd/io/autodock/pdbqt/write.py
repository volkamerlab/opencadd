from typing import Sequence, Tuple, Union, Optional
from pathlib import Path
import numpy as np
import pandas as pd
from openbabel import pybel
import opencadd as oc
from opencadd import _typing


def from_ensemble(
        ensemble,
        output_filename: Optional[str] = None,
        output_path: _typing.PathLike = None,
        models: Optional[Union[int, Sequence[int]]] = None,
):
    pdb_strings = oc.io.pdb.write.from_chemsys(
        system=ensemble,
        models=models,
        separate_models=True
    )
    return_vals = []
    temp_filepath = Path.cwd() / "_temp_opencadd.pdb"
    for pdb_string in pdb_strings:
        with open(temp_filepath, "wt") as f:
            f.write(pdb_string)
        return_vals.append(
            from_pdb_filepath(filepath=temp_filepath)
        )
    return return_vals


def from_pdb_filepath(
        filepath: _typing.PathLike,
        output_filename: Optional[str] = None,
        output_path: Optional[_typing.PathLike] = None,
):
    """
    Convert a PDB file to a PDBQT file, and save it in the given filepath.

    Parameters
    ----------
    filepath: str or pathlib.Path
        Path to input PDB file.
    output_path: str or pathlib.Path
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

    References
    ----------
    https://open-babel.readthedocs.io/en/latest/FileFormats/AutoDock_PDBQT_format.html
    """
    # pybel.readfile() provides an iterator over the Molecules in a file.
    # To access the first (and possibly only) molecule in a file, we use next()
    input_path = Path(filepath)
    molecule = next(pybel.readfile("pdb", str(input_path)))
    # if protonate_for_pH:
    molecule.OBMol.CorrectForPH(7.4)
    molecule.addh()
    # if add_hydrogens:

    # if calculate_partial_charges:
    for atom in molecule.atoms:
        atom.OBAtom.GetPartialCharge()
    # TODO: expose write options to function sig (see ref.)
    if output_path is None:
        output_filepath = (input_path.parent / input_path.stem).with_suffix(".pdbqt")
    else:
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        output_filepath = (output_path/output_filename).with_suffix(".pdbqt")
    molecule.write(
        format="pdbqt",
        filename=str(output_filepath),
        overwrite=True,
        opt={"r": None, "n": None, "p": None, "h":None}
    )
    if output_filename is not None:
        return output_filepath.resolve()
    with open(output_filepath) as f:
        pdbqt_str = f.read()
    output_filepath.unlink()
    return pdbqt_str
