"""
Functions for handling pdbqt files.
"""

# Standard library
import os
from typing import NoReturn, Tuple, Optional, Union
from pathlib import Path
# 3rd-party
import numpy as np
import pandas as pd
from openbabel import pybel
# Self
from opencadd.misc.parsing import extract_column_from_string_array
from opencadd.consts import autodock


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


def extract_autodock_atom_types_from_pdbqt(filepath_pdbqt: Path) -> Tuple[str]:
    """
    Extract the AutoDock-defined atom types for all atoms in a PDBQT file.

    Parameters
    ----------
    filepath_pdbqt : pathlib.Path
        Path to the PDBQT file.

    Returns
    -------
    tuple[str]
        Type of each atom in the same order as in the file, where each type is defined by AutoDock
        as a string of one or two characters, such as "A", "C", "HD", "N", "NA", "OA", "SA", etc.
    """
    # TODO: Create a complete PDBQT parser class and add this and other related functions as
    #  methods.
    with open(filepath_pdbqt, "r") as f:
        lines = f.readlines()
    atom_types = [line.split()[-1] for line in lines if line.startswith("ATOM")]
    return tuple(set(atom_types))


def parse_pdbqt(filepath_pdbqt: Path):
    """
    Parse a PDBQT file.

    Parameters
    ----------
    filepath_pdbqt : pathlib.Path
        Path to the PDBQT file.

    Returns
    -------
    dict(str, pandas.DataFrame)
        A dictionary of record names (e.g. "ATOM") and their corresponding dataframes.

    References
    ----------
    PDB file format documentation:
        https://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf
    """

    def parse_atom_records(record_lines_atom: np.ndarray):
        """
        Parse ATOM records
        """
        columns = {
            "serial": ((7, 11), int),
            "name": ((13, 16), (str, 4)),
            "altLoc": ((17, 17), (str, 1)),
            "resName": ((18, 20), (str, 3)),
            "chainID": ((22, 22), (str, 1)),
            "resSeq": ((23, 26), int),
            "iCode": ((27, 27), (str, 1)),
            "x": ((31, 38), float),
            "y": ((39, 46), float),
            "z": ((47, 54), float),
            "occupancy": ((55, 60), float),
            "tempFactor": ((61, 66), float),
            "partial_charge": ((67, 76), float),
            "autodock_atom_type": ((78, 79), (str, 2))
        }

        df = pd.DataFrame()
        for col_name, (col_range, col_dtype) in columns.items():
            df[col_name] = np.char.strip(
                extract_column_from_string_array(
                    array=record_lines_atom,
                    char_range=(col_range[0] - 1, col_range[1])
                )
            ).astype(col_dtype)

        autodock_atom_types = autodock.AtomTypes()
        autodock_atom_types_data = [
            np.array([getattr(atom_type, attr) for atom_type in autodock_atom_types])
            for attr in ["name", "hbond_status", "hbond_count"]
        ]
        indices_target_atom_types = np.where(
            df.autodock_atom_type.values[..., np.newaxis] == autodock_atom_types_data[0]
        )[1]
        df["hbond_status"] = autodock_atom_types_data[1][indices_target_atom_types]
        df["hbond_count"] = autodock_atom_types_data[2][indices_target_atom_types]
        return df

    record_parsers = {
        "ATOM" : parse_atom_records
    }

    with open(filepath_pdbqt, "r") as f:
        lines = np.array(f.readlines())

    records = dict()
    for record, parser in record_parsers.items():
        record_mask = np.char.startswith(lines, prefix=record)
        records[record] = parser(lines[record_mask])
    return records
