"""
T2F-Pharm (Truly Target-Focused Pharmacophore) model.
"""


# Standard library
from typing import Sequence, Union, Optional, Literal, Tuple, Dict
from pathlib import Path
# 3rd-party
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
# Self
from opencadd.docking import autogrid
from opencadd.io import pdbqt
from opencadd.io.pdbqt import TYPE_AUTODOCK_ATOM_TYPE, DATA_AUTODOCK_ATOM_TYPE
from opencadd.bindingsite.utils import count_psp_events_from_vacancy_tensor
from opencadd.misc.spatial import Grid


class T2FPharm:
    """
    Truly Target Focused (T2F) pharmacophore model.
    """

    def __init__(
            self,
            filepath_target_structure: Union[str, Path],
            path_output: Optional[Union[str, Path]] = None,
            type_probes: Sequence[TYPE_AUTODOCK_ATOM_TYPE] = ("A", "C", "HD", "OA"),
            coords_grid_center: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
            spacing_grid_points: float = 0.6,
            dims_grid: Optional[Tuple[float, float, float]] = (16, 16, 16),
            npts: Optional[Tuple[int, int, int]] = None,
    ):
        """
        Get interaction energies between the target structure and each of the input probe types
        (plus electrostatic and desolvation energies), calculated using AutoGrid4 at every point
        on a regular (i.e. evenly spaced) cuboid grid, placed on (usually the binding pocket of)
        the target structure, as defined by the input arguments.


        Parameters
        ----------
        filepath_target_structure : str | pathlib.Path
            Path to the input PDBQT file containing the target structure.
        path_output : str | pathlib.Path, Optional, default: None
            An output path to store the results. If not provided, the output files will be
            stored in the same folder as the input file. If a non-existing path is given,
            a new directory will be created with all necessary parent directories.
        type_probes : Sequence[Literal]
            AutoDock-defined atom types (e.g. "A", "C", "HD", "N", "NA", "OA", "SA", "Cl") of
            the probes, for which interaction energies must be calculated.
            For a full list of atom types, see `opencadd.io.pdbqt.DATA_AUTODOCK_ATOM_TYPE`.
        coords_grid_center : tuple[float, float, float] | "auto"
            Coordinates (x, y, z) of the center of grid, in the reference frame of the target
            structure, in Ångstrom (Å). If set to "auto", AutoGrid automatically centers
            the grid on the target structure's center.
        spacing_grid_points : float, Optional, default: 0.6
            The grid-point spacing, i.e. distance between two adjacent grid points in Å.
            Grid points are orthogonal and uniformly spaced in AutoGrid4, i.e. this value is
            used for all three dimensions.
        dims_grid : tuple[float, float, float], Optional, default: (16, 16, 16)
            Length of the grid (in Å), along x-, y-, and z-axis, respectively.
        npts : tuple[float, float, float], Optional, default: None
            Number of grid points along x-, y-, and z-axis respectively (loosely defined; for an
            accurate definition, see function `opencadd.docking.autogrid.routine_run_autogrid`).
            If not provided, this will be automatically calculated from input arguments
            `dims_grid` and `spacing_grid_points`, since AutoGrid4 actually requires
            `coords_grid_center`, `spacing_grid_points` and `npts`, in order to define a grid.
            Notice that since `npts` values can only be even integers, calculating `npts` may
            result in a grid that is slightly larger than the input dimensions (for more info,
            see the function `opencadd.docking.autogrid.calculate_npts`). Therefore, if this
            argument is provided, the `dims_grid` argument will be ignored.
        """
        # Verify the input file is an existing PDBQT file and assign output path.
        filepath_target = Path(filepath_target_structure)
        if not filepath_target.is_file():
            raise ValueError("Input filepath does not point to a file.")
        try:
            df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=filepath_target)["ATOM"]
        except Exception as e:
            raise ValueError("Could not parse the input PDBQT file.") from e
        if path_output is None:
            path_output = filepath_target.parent
        else:
            path_output = Path(path_output)
            path_output.mkdir(parents=True, exist_ok=True)
        # Verify that all input probe types are valid AutoDock atom types
        type_probes = np.array(type_probes)
        probe_is_invalid = np.isin(
            element=type_probes,
            test_elements=DATA_AUTODOCK_ATOM_TYPE.type.values,
            invert=True
        )
        if np.any(probe_is_invalid):
            raise ValueError(
                f"The following `type_probes` are not valid AutoDock atom types: "
                f"{type_probes[probe_is_invalid]}"
            )
        # Set instance attributes.
        self._filepath_target: Path = filepath_target
        self._path_output: Path = path_output
        self._data_target_structure: pd.DataFrame = df_pdbqt
        # Calculate grid energies using AutoGrid4
        self._grid: Grid = autogrid.routine_run_autogrid(
            receptor=self._filepath_target,
            gridcenter=coords_grid_center,
            npts=npts if npts is not None else autogrid.calculate_npts(
                dimensions_pocket=dims_grid,
                spacing=spacing_grid_points
            ),
            spacing=spacing_grid_points,
            ligand_types=type_probes,
            path_output=self._path_output
        )
        return

    @property
    def grid(self) -> Grid:
        """
        The complete grid.

        Returns
        -------
        numpy.ndarray
        A 4-dimensional array, where the first three dimensions represent the grid points,
        and the last dimension contains the data for each grid point. It is the direct
        return value of the function `opencadd.docking.autogrid.routine_run_autogrid`;
        for more information, see the documentation of that function.
        """
        return self._grid

    @property
    def target(self) -> pd.DataFrame:
        return self._data_target_structure

    def calculate_mask_vacant(
            self,
            energy_refs: Sequence[Union[TYPE_AUTODOCK_ATOM_TYPE, Literal["e", "d"]]],
            energy_cutoff: float = +0.6,
            mode: Optional[Literal["max", "min", "avg", "sum"]] = "sum",
    ) -> np.ndarray:
        """
        Calculate whether each grid point is vacant (or occupied by a target atom), based on given
        reference energy types, and a cutoff value.

        Parameters
        ----------
        energy_refs : Sequence[Literal]
            Grid-point energy value(s) to take as reference. These can be from the probe
            types for which energy grids were calculated (i.e. elements of input argument
            `type_probes` in method `calculate_energy_grid`), plus "e" and "d", referring to the
            electrostatic potential, and desolvation energy grids, respectively.
        energy_cutoff : float, Optional, default: +0.6
            Cutoff value for energy; grid points with energies lower than cutoff are considered
            vacant.
        mode: Literal["max", "min", "avg", "sum"], Optional, default: "sum"
            If more than one energy type is inputted in `energy_ref`, this parameter defines how
            those different energy values must be processed, before comparing with the cutoff
            value. If only one energy type is entered, then this parameter is ignored.

        Returns
        -------
        mask_vacant : numpy.ndarray
            A 3-dimensional boolean array matching the dimensions of the energy grid,
            indicating whether each grid point is vacant (True), or occupied (False).
            Thus, the vacant grid points can easily be indexed using boolean indexing with this
            array, i.e.: `grid[mask_vacant]`.
        """
        # The reducing operations corresponding to each `mode`:
        reducing_op = {"max": np.max, "min": np.min, "avg": np.mean, "sum": np.sum}
        # Verify that all given reference energy types are already calculated
        not_calculated = np.isin(energy_refs, self._type_energy_grids, invert=True)
        if np.any(not_calculated):
            raise ValueError(
                f"The following energy grids were not calculated: {energy_refs[not_calculated]}"
            )
        # Take the subgrid corresponding to given energy references
        subgrid_ref = self._grid.get_properties(names=energy_refs)
        # If only one energy reference is selected, take that, otherwise reduce the given
        # references using the given operation.
        energy_vals = subgrid_ref if len(energy_refs) == 1 else reducing_op[mode](subgrid_ref, axis=-1)
        # Apply cutoff and return
        return energy_vals < energy_cutoff

    def calculate_mask_buried(
            self,
            mask_vacancy: np.ndarray,
            num_directions: Literal[3, 7, 13] = 7,
            psp_len_max: float = 10.0,
            psp_count_min: int = 4,
    ) -> np.ndarray:
        """
        Calculate whether each grid point is buried inside the target structure or not, based on
        counting the number of protein-solvent-protein (PSP) events for each point, and applying
        a cutoff.

        Notice that before using this method, a vacancy grid must be calculated using the method
        `calculate_vacancy_from_energy`.

        Parameters
        ----------
        num_directions : Literal[3, 7, 13], Optional, default: 7
            Number of directions to look for PSP events for each grid point. Directions are
            all symmetric around each point, and each direction entails both positive and
            negative directions, along an axis. Directions are defined in a 3x3x3 unit cell,
            i.e. for each point, there are 26 neighbors, and thus max. 13 directions.
            Options are:
            3: x-, y- and z-directions
            7: x, y, z, and four 3d-diagonals (i.e. with absolute unit vector [1, 1, 1])
            13: x, y, z, four 3d-diagonals, and six 2d-diagonals (e.g. [1, 1, 0])
        psp_len_max : float, Optional, default: 10.0
            Maximum acceptable distance for a PSP event, in Ångstrom (Å).
        psp_count_min : int, Optional, default: 4
            Minimum required number of PSP events for a grid point, in order to count as buried.

        Returns
        -------
        mask_buriedness : numpy.ndarray, Optional, default: None
            A 3-dimensional array matching the dimensions of the energy grid,
            indicating whether each grid point is buried (True), or exposed (False).
            Thus, the buried grid points can easily be indexed using boolean indexing with this
            array, i.e.: `grid[mask_buried]`.
        """
        grid_psp_counts = count_psp_events_from_vacancy_tensor(
            grid_vacancy=mask_vacancy,
            spacing_grid_points=self._grid.spacing,
            num_directions=num_directions,
            psp_len_max=psp_len_max
        )
        mask_buried = np.zeros(shape=self._grid.shape, dtype=np.bool_)
        mask_buried[mask_vacancy] = grid_psp_counts >= psp_count_min
        return mask_buried

    def count_hbond_partners_in_target(
            self,
            max_len_hbond: float = 3.0,
            mask_grid: Optional[np.ndarray] = None,
    ):
        if mask_grid is None:
            mask_grid = np.ones(shape=self._grid.shape, dtype=np.bool_)
        coords_target_atoms = self._data_target_structure[["x", "y", "z"]].to_numpy()
        mask_proximate = self._grid.calculate_mask_proximate(
            coords_targets=coords_target_atoms,
            distance=max_len_hbond,
            mask=mask_grid,
        )
        count_acc_don_in_target = np.zeros(shape=(*self._grid.shape, 2), dtype=np.byte)
        for idx, hbond_role in enumerate(["is_acceptor", "is_donor"]):
            count_acc_don_in_target[mask_grid, idx] = np.count_nonzero(
                np.logical_and(
                    self._data_target_structure[hbond_role],
                    mask_proximate
                ),
                axis=-1
            )
        return count_acc_don_in_target

    def calculate_mask_feature(
            self,
            mask_grid: np.ndarray,
            count_acceptors_in_target: np.ndarray,
            max_energy: float = -0.35,
            min_count_acceptor_atoms: int = 1
    ):
        mask_feature_hbond_donor = np.zeros(shape=self._grid.shape[:3], dtype=np.bool_)
        mask_feature_hbond_donor[mask_grid] = np.logical_and(
            self._grid[mask_grid, 2] <= max_energy,
            count_acceptors_in_target[mask_grid] >= min_count_acceptor_atoms
        )
        return mask_feature_hbond_donor


    def calculate_features(
            self,
            mask_vacant: np.ndarray,
            mask_buried: np.ndarray,
            mask_proximate: np.ndarray,
            count_hbond_acc: np.ndarray,
            count_hbond_donor: np.ndarray,
            cutoff_energies: Tuple[
                Tuple[Union[TYPE_AUTODOCK_ATOM_TYPE, Literal["e", "d"]], float]
            ] = (("A", -0.4), ("C", -0.4), ("HD", -0.35), ("OA", -0.6), ("e", -1.2)),
    ):



# def routine_run_t2fpharm(
#             filepath_target_structure: Union[str, Path],
#             path_output: Optional[Union[str, Path]] = None,
#             type_probes: Sequence[TYPE_AUTODOCK_ATOM_TYPE] = ("A", "C", "HD", "OA"),
#             coords_grid_center: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
#             spacing_grid_points: float = 0.6,
#             dims_grid: Optional[Tuple[float, float, float]] = (16, 16, 16),
#             npts: Optional[Tuple[int, int, int]] = None,
#             energy_refs: Sequence[Union[TYPE_AUTODOCK_ATOM_TYPE, Literal["e", "d"]]] = ["OA"],
#             energy_cutoff: float = +0.6,
#             mode: Optional[Literal["max", "min", "avg", "sum"]] = "sum",
# ):
#     t2fpharm = T2FPharm(
#         filepath_target_structure=filepath_target_structure,
#         path_output=path_output,
#         type_probes=type_probes,
#         coords_grid_center=coords_grid_center,
#         spacing_grid_points=spacing_grid_points,
#         dims_grid=dims_grid,
#         npts=npts,
#     )
#
#     mask_vacant = t2fpharm.calculate_mask_vacant(
#         energy_refs=energy_refs,
#         energy_cutoff=energy_cutoff,
#         mode=mode
#     )
#
#     mask_buried = t2fpharm.calculate_mask_buried(
#         mask_vacancy=mask_vacant,
#         num_directions=num_directions,
#         psp_len_max=,
#         psp_count_min=
#     )
#
#     mask_proximate = t2fpharm.calculate_mask_proximate(
#         mask_grid=
#     )
#
#
#
#
#
#
#         is_donor = DATA_AUTODOCK_ATOM_TYPE.is_donor.values[indices_target_atom_types]
#         is_acceptor = DATA_AUTODOCK_ATOM_TYPE.is_acceptor.values[indices_target_atom_types]
#
#
#
#
#         atom_types_hbond_donor = np.array(["HD", "HS"])
#         atom_types_hbond_acceptor = np.array(["NA", "NS", "OA", "OS", "SA"])
#         grid_acceptor_donor_count = []
#         for grid_point_proximity_mask in proximity_mask:
#             proximate_atom_types = pdbqt_df.autock_atom_type.values[grid_point_proximity_mask]
#             acceptor_count = np.count_nonzero(
#                 np.isin(proximate_atom_types, atom_types_hbond_acceptor))
#             donor_count = np.count_nonzero(np.isin(proximate_atom_types, atom_types_hbond_donor))
#             grid_acceptor_donor_count.append([acceptor_count, donor_count])
#         grid_acceptor_donor_count = np.array(
#             grid_acceptor_donor_count, dtype=np.int8
#         ).reshape((*subgrid.shape[:3], 2))
#         hbond_mask = grid_acceptor_donor_count.sum(axis=-1) > 0
#
#         points_hydrophilic = subgrid[hbond_mask]
#         points_hydrophobic = subgrid[~hbond_mask]
#
#         points_hbond_donor_mask = (grid_acceptor_donor_count[..., 0] > 0) & (subgrid)


# def t2fpharm(
#         filepath_target_structure: Path,
#         coordinates_pocket_center: Union[Tuple[float, float, float], str] = "auto",
#         dimensions_pocket: Tuple[float, float, float] = (16, 16, 16),
#         spacing_grid: float = 0.6,
#         atom_types_target: Sequence[str] = ("A", "C", "HD", "N", "NA", "OA", "SA", "Cl"),
#         types_maxpot_probes: Sequence[Tuple[str, float]] = (
#                 ("A", -0.4),
#                 ("C", -0.4),
#                 ("HD", -0.35),
#                 ("OA", -0.6),
#         ),
#         maxpot_abs_electrostat: float = 1.0,
#         maxpot_occupancy: float = +0.6,
#         count_dirs_psp: Literal[3, 7, 13] = 7,
#         max_psp_len: float = 10.0,
#         min_psp_count: int = 4,
#         max_hbond_len: float = 3.0,
#         min_neighbor_dist_clustering: float = 1.21,
#         min_common_neighbor_count_clustering: int = 6,
#         min_points_per_cluster_count: int = 15,
#         path_output: Optional[Path] = None
# ):

a=np.array([1,2,3])
print(a[1])

pharm = T2FPharm(filepath_target_structure="/Users/home/Downloads/test_run/3w32.pdbqt")
mask_vacant = pharm.calculate_mask_vacant(
    energy_refs=["e","d"],
)
pharm.calculate_mask_buried(mask_vacant)