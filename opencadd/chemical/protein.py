from typing import Sequence, Tuple, Union, Optional
from pathlib import Path
import numpy as np
import pandas as pd
import scipy as sp
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
        self._trajectory_2d_view: np.ndarray = self._trajectory.reshape(-1, 3)
        self._kdtrees: Tuple[sp.spatial.KDTree] = None
        self._kdtree: sp.spatial.KDTree = sp.spatial.KDTree(data=self._trajectory_2d_view)
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
    def kdtree(self):
        return self._kdtree

    @property
    def kdtrees(self):
        if self._kdtrees is None:
            self._kdtrees = tuple(
                sp.spatial.KDTree(data=self._trajectory[i]) for i in range(self.trajectory_length)
            )
        return self._kdtrees

    def nearest_atoms(
            self,
            coords: np.ndarray,
            num_atoms: int = 1,
            error_tolerance: float = 0,
            trajectory_filter: slice = slice(None),
    ):
        """
        For a given number of points, find the nearest atom in the protein to each point.

        Parameters
        ----------
        coords : ndarray, shape: (n, 3)
            Cartesian coordinates of n points, for which the nearest atom in protein must be found.
        num_atoms : int, optional, default: 1
            Number of the nearest atoms to find for each point.
        error_tolerance :
        trajectory_filter : slice, optional, default: slice(None)
            Slice of trajectories to consider.

        Returns
        -------
        distances, indices : ndarray, ndarray
            Distances to, and indices of the k nearest atoms to each point, in each trajectory.
            For t selected trajectories, the shape of each array will be (t, n) when `num_atoms`
            is 1, or (t, n, `num_atoms`) otherwise, with n being the number of points in `coords`.
        """
        num_selected_trajectories = len(self.kdtrees[trajectory_filter])
        distances = np.empty(
            shape=(num_selected_trajectories, coords.shape[0] * num_atoms),
            dtype=np.single
        )
        indices = np.empty(
            shape=(num_selected_trajectories, coords.shape[0] * num_atoms),
            dtype=np.ushort
        )
        for traj_idx, kdtree in enumerate(self._kdtrees[trajectory_filter]):
            nearest_atom_distances, nearest_atom_indices = kdtree.query(
                coords,
                k=num_atoms,
                eps=error_tolerance,
                workers=-1,
            )
            distances[traj_idx] = nearest_atom_distances
            indices[traj_idx] = nearest_atom_indices
        return distances, indices

    def distance_to_atoms(
            self,

    ):
        pass

    def distance_to_atoms_sparse(
            self,
            kdtree: sp.spatial.KDTree,
            max_distance: float,
    ):
        distances = self._kdtree.sparse_distance_matrix(
            other=kdtree,
            max_distance=max_distance
        ).toarray().reshape(self.trajectory_length, self.count_atoms, -1)
        return np.moveaxis(distances, source=2, destination=1)

    @property
    def count_atoms(self) -> int:
        return len(self._atom_data)

    @property
    def atom_data(self):
        return self._atom_data

    @property
    def autodock_atom_types(self) -> np.ndarray:
        return

    @property
    def trajectory_length(self) -> int:
        return self._trajectory.shape[0]

    @property
    def trajectory(self) -> np.ndarray:
        return self._trajectory

    @property
    def view(self):
        return nglview.show_file(str(self._structure_filepath))



    def create_new_ngl_widget(self):
        return nglview.show_file(str(self._structure_filepath), height="800px")


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


def from_pdb_id(pdb_id: str) -> Protein:
    pass


def from_filepath(path: PathLike) -> Protein:
    pass


def from_file_content(content: Union[str, bytes]) -> Protein:
    pass


def to_file_pdbqt(protein: Protein) -> Path:
    pass


