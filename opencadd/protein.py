from typing import Sequence, Tuple, Union, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import nglview
from openbabel import pybel

from opencadd.typing import PathLike
from opencadd.misc.parsing import extract_column_from_string_array
from opencadd.consts import autodock


class Protein:

    _SUPPORTED_INPUT_FILES = ("pdb", "pdbqt")

    def __init__(
            self,
            name: str = None,
            pdb_code: str = None,
            atom_data: pd.DataFrame = None,
            trajectory: np.ndarray = None,
    ):

        self._name = name
        self._pdb_code = pdb_code
        self._atom_data = atom_data
        self._trajectory = trajectory
        return

    @classmethod
    def from_file(cls, paths: Union[PathLike, Sequence[PathLike]]):
        """
        Instantiate a Protein from file.

        Parameters
        ----------
        path : PathLike
            Path

        Returns
        -------

        """
        if isinstance(paths, PathLike):
            paths = [paths]
        filepaths = [Path(path) for path in paths]
        filetype = filepaths[0].suffix.lower()[1:]
        if filetype not in Protein._SUPPORTED_INPUT_FILES:
            raise ValueError(f"File type '{filetype}' is not supported.")
        return getattr(Protein, f"_from_file_{filetype}")(paths=filepaths)

    @property
    def atom_data(self):
        return self._atom_data

    @property
    def autodock_atom_types(self) -> np.ndarray:
        return

    @property
    def trajectory(self) -> np.ndarray:
        return self._trajectory

    @property
    def view(self):
        return nglview.show_file(str(self._structure_filepath))

    def to_file_pdbqt(self):
        pass


    @classmethod
    def _from_file_pdb(cls, paths: Sequence[Path]):
        pass

    @classmethod
    def _from_file_pdbqt(cls, paths: Sequence[Path]):
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

        with open(paths[0], "r") as f:
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

    def pdb_to_pdbqt_openbabel(
            self,
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

        References
        ----------
        https://open-babel.readthedocs.io/en/latest/FileFormats/AutoDock_PDBQT_format.html
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
        # TODO: expose write options to function sig (see ref.)
        molecule.write(
            format="pdbqt",
            filename=str(pdbqt_filepath.with_suffix(".pdbqt")),
            overwrite=True,
            opt={"r": None, "n": None, "p": None}
        )
        return molecule
