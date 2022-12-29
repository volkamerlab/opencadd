from typing import Sequence
import numpy as np
import pandas as pd

from opencadd.typing import PathLike


def read_file(filepath: PathLike):
    """
    Parse a PDBQT file.

    Parameters
    ----------
    filepath : PathLike
        Path to PDB file.

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

        autodock_atom_types_ids = np.array(
            [atom_type.name for atom_type in autodock.AtomType])
        autodock_atom_types_data = [
            np.array([getattr(atom_type, attr) for atom_type in autodock.AtomType])
            for attr in ["hbond_status", "hbond_count"]
        ]
        indices_target_atom_types = np.where(
            df.autodock_atom_type.values[..., np.newaxis] == autodock_atom_types_ids
        )[1]
        df["hbond_acc"] = autodock_atom_types_data[0][indices_target_atom_types] == 1
        df["hbond_don"] = autodock_atom_types_data[0][indices_target_atom_types] == -1
        df["hbond_count"] = autodock_atom_types_data[1][indices_target_atom_types]
        return df

    record_parsers = {
        "ATOM": parse_atom_records
    }

    with open(filepath[0], "r") as f:
        lines = np.array(f.readlines())

    records = dict()
    for record, parser in record_parsers.items():
        record_mask = np.char.startswith(lines, prefix=record)
        records[record] = parser(lines[record_mask])

    trajectory = records["ATOM"][["x", "y", "z"]].to_numpy()[np.newaxis]

    return cls(
        atom_data=records["ATOM"],
        trajectory=trajectory,
    )
