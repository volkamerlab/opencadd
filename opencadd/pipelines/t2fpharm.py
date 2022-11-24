"""
T2F-Pharm (Truly Target-Focused Pharmacophore) model.
"""


# Standard library
from typing import Sequence, Union, Optional, Literal, Tuple, Dict
from pathlib import Path
# 3rd-party
import numpy as np
import pandas as pd
# Self
from opencadd.docking import autogrid
from opencadd.io import pdbqt
from opencadd.consts import autodock
from opencadd.misc.spatial import grid_distance, Grid

class T2FPharm(Grid):
    """
    Truly Target Focused (T2F) pharmacophore model.

    A 3-dimensional array containing the calculated energy values for each grid point.
        The shape of the array is (nx, ny, nz), where nx, ny, and nz are the number of grid
        points along the x-, y-, and z-axis, respectively. These are each equal to the
        corresponding input `npts` value, plus 1 (due to center point). The array is ordered in
        the same way that the grid points are ordered in space, and thus the individual grid
        points can be indexed using their actual coordinates (in unit vectors), assuming the
        origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
        grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.

    A 4-dimensional array, where the first three dimensions represent the grid points,
        and the last dimension contains the data for each grid point.

        Each grid point contains energy values for each of the input ligand atom-types (probes),
        in the given order. Thus, with 'nt' being the number of input probes, the first 'nt'
        elements of the last dimension of the grid array correspond to their energy values.
        There are five other elements at the end of the last dimension, namely the electrostatic
        potential and desolvation energy, followed by the x-, y- and z-coordinates of the grid
        point in the reference frame of the target structure.

        The shape of the array is thus (nx, ny, nz, nt + 5), where nx, ny, and nz are the number of
        grid points along the x-, y-, and z-axis, respectively. These are each equal to their
        corresponding `npts` value, plus 1 (due to center point). The array is ordered in
        the same way that the grid points are ordered in space, and thus the individual grid
        points can be indexed using their actual coordinates (in unit vectors), assuming the
        origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
        grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.

    Get interaction energies between the target structure and each of the input probe types
        (plus electrostatic and desolvation energies), calculated using AutoGrid4 at every point
        on a regular (i.e. evenly spaced) cuboid grid, placed on (usually the binding pocket of)
        the target structure, as defined by the input arguments.
    """

    def __init__(
            self,
            receptor_filepath: Path,
            ligand_types: Sequence[autodock.AtomType],
            grid_center: Union[Tuple[float, float, float], Literal["auto"]],
            grid_size: Tuple[float, float, float] = (16, 16, 16),
            grid_spacing: float = 0.6,
            smooth: float = 0.5,
            dielectric: float = -0.1465,
            param_filepath: Optional[Path] = None,
            output_path: Optional[Path] = None,
    ):
        """

        Parameters
        ----------
        receptor_filepath : pathlib.Path
            Path to the PDBQT structure file of the macromolecule.
        ligand_types : Sequence[opencadd.consts.autodock.AtomType]
            Types of ligand atoms, for which interaction energies must be calculated.
            For more information, see `opencadd.consts.autodock.AtomType`.
        grid_center : tuple[float, float, float] | Literal["auto"], Optional, default: "auto"
            Coordinates (x, y, z) of the center of grid map, in the reference frame of the target
            structure, in Ångstrom (Å). If set to "auto", AutoGrid automatically centers the grid
            on the receptor's center of mass.
        grid_npts : Tuple[int, int, int], Optional, default: (40, 40, 40)
            Number of grid points to add to the central grid point, along x-, y- and z-axes,
            respectively. Each value must be an even integer number; when added to the central grid
            point, there will be an odd number of points in each dimension. The number of x-, y and
            z-grid points need not be equal.
        grid_spacing : float, Optional, default: 0.375
            The grid-point spacing, i.e. distance between two adjacent grid points in Ångstrom (Å).
            Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
            all three dimensions.
        smooth : float, Optional, default: 0.5
            Smoothing parameter for the pairwise atomic affinity potentials (both van der Waals
            and hydrogen bonds). For AutoDock4, the force field has been optimized for a value of
            0.5 Å.
        dielectric : float, Optional, default: -0.1465
            Dielectric function flag: if negative, AutoGrid will use distance-dependent dielectric
            of Mehler and Solmajer; if the float is positive, AutoGrid will use this value as the
            dielectric constant. AutoDock4 has been calibrated to use a value of –0.1465.
        param_filepath : pathlib.Path, Optional, default: None
            User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.
        output_path: pathlib.Path
            Path to a folder to write the output files in. If not provided, the output files will be
            stored in the same folder as the input file. If a non-existing path is given,
            a new directory will be created with all necessary parent directories.
        """
        npts = np.array(autogrid.calculate_npts(grid_size=grid_size, grid_spacing=grid_spacing))
        energies, filepaths_mapfiles = autogrid.routine_run(
            receptor_filepath=receptor_filepath,
            ligand_types=ligand_types,
            grid_center=grid_center,
            grid_npts=npts,
            grid_spacing=grid_spacing,
            smooth=smooth,
            dielectric=dielectric,
            param_filepath=param_filepath,
            output_path=output_path,
        )
        grid_shape = npts + 1
        if grid_center == "auto":
            grid_center, _, _ = autogrid.extract_grid_params(map_filepath=filepaths_mapfiles[0])
        super().__init__(
            shape=grid_shape,
            coords_origin=np.array(grid_center) - grid_spacing * npts / 2,
            spacing=grid_spacing,
            data=[energy_value.reshape(tuple(grid_shape), order="F") for energy_value in energies],
        )
        self._ligand_types = np.array(ligand_types, dtype=object)
        self._path_mapfiles = filepaths_mapfiles
        self._output_path = receptor_filepath.parent if output_path is None else output_path
        self._output_path.mkdir(parents=True, exist_ok=True)
        self._df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptor_filepath)["ATOM"]
        return

    @property
    def receptor_data(self) -> pd.DataFrame:
        return self._df_pdbqt

    @property
    def grid_electrostat_pot(self) -> np.ndarray:
        """
        Grid of calculated electrostatic potentials.
        """
        return self._tensor[..., -2]

    @property
    def grid_desolvation(self) -> np.ndarray:
        """
        Grid of calculated desolvation energies.
        """
        return self._tensor[..., -1]

    @property
    def grid_probes(self) -> np.ndarray:
        """
        Grid of calculated probe interaction energies.
        """
        return self._tensor[..., :-2]

    def grid_probe(self, ligand_type: autodock.AtomType) -> np.ndarray:
        """
        Grid of calculated probe interaction energies.
        """
        idx = np.argwhere(self._ligand_types == ligand_type)
        if idx.size == 0:
            raise IndexError("No energy has been calculated for the given ligand type.")
        return self._tensor[..., idx[0, 0]]

    def calculate_mask_vacant(
            self,
            energy_refs: Sequence[Literal["e", "d"]],
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
        grid_psp_counts = self.count_psp_events(
            mask_vacancy=mask_vacancy,

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
        coords_target_atoms = self._df_pdbqt[["x", "y", "z"]].to_numpy()
        mask_proximate = self.calculate_mask_proximate(
            coords_targets=coords_target_atoms,
            distance=max_len_hbond,
            mask=mask_grid,
        )
        counts = []
        for hbond_role in ["hbond_acc", "hbond_don"]:
            count = np.zeros(shape=self.shape, dtype=np.byte)
            count[... if mask_grid is None else mask_grid] = np.count_nonzero(
                np.logical_and(
                    self._df_pdbqt[hbond_role],
                    mask_proximate
                ),
                axis=-1
            )
            counts.append(count)
        return tuple(counts)

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


    def count_psp_events(
            self,
            mask_vacancy: np.ndarray,
            num_directions: Literal[3, 7, 13] = 7,
            psp_len_max: float = 10.0,
    ) -> np.ndarray:
        """
        Count the number of protein-solvent-protein (PSP) events for each point on a 3-dimensional
        grid, given a 3-dimensional boolean array indicating the vacancy of each grid point.

        Parameters
        ----------
        mask_vacancy : numpy.ndarray
            A 3-dimensional boolean array representing a 3-dimensional grid, indicating
            whether each grid point is vacant (True), or occupied (False).
        num_directions : Literal[3, 7, 13]
            Number of directions to look for PSP events for each grid point. Directions are
            all symmetric around each point, and each direction entails both positive and
            negative directions, along an axis. Directions are defined in a 3x3x3 unit cell,
            i.e. for each point, there are 26 neighbors, and thus max. 13 directions.
            Options are:
            3: x-, y- and z-directions
            7: x, y, z, and four 3d-diagonals (i.e. with absolute direction unit vector [1, 1, 1])
            13: x, y, z, four 3d-diagonals, and six 2d-diagonals (e.g. [1, 1, 0])
        psp_len_max : float
            Maximum acceptable distance for a PSP event, in the same units as `spacing_grid_points`.

        Returns
        -------
        grid_psp_count : numpy.ndarray
            A 3-dimensional array with the same shape as the input `grid_vacancy` argument,
            where each value corresponds to the number of found PSP events for that grid point.
        """
        if num_directions not in [3, 7, 13]:
            raise ValueError("Argument `num_directions` must be either 3, 7, or 13.")
        # Get all 26 half directions, in the order: orthogonal, 3d-diagonal, 2d-diagonal
        dirs = np.concatenate([self._DIRS_ORTHO, self._DIRS_DIAG_3D, self._DIRS_DIAG_2D])
        # Calculate the length of each direction unit vector
        len_dir_vectors = np.repeat(self.spacing, 26) * np.repeat(
            np.array([1, np.sqrt(3), np.sqrt(2)]),
            repeats=[6, 8, 12],
        )
        # Calculate the maximum allowed number of times to travel in each half-direction
        max_mult_factors = psp_len_max // len_dir_vectors
        # Calculate distance of each vacant grid point to the nearest occupied grid point in each
        # half direction, in units of corresponding distance vectors
        dists_single_dir = grid_distance(
            grid=mask_vacancy,
            start_indices=np.argwhere(mask_vacancy),
            directions=dirs[:num_directions * 2],
            target_value=0,
            max_dist_in_dir=max_mult_factors[:num_directions * 2]
        )
        # For each grid-point, filter directions where one half-direction is zero (meaning no
        # neighbor was found in that direction).
        has_psp = np.logical_and(dists_single_dir[..., 0::2] != 0,
                                 dists_single_dir[..., 1::2] != 0)
        # Add distances to neighbors in positive half-directions , to distances to neighbors in
        # negative half-directions, in order to get the PSP length in units of direction vectors,
        # and then multiply by direction unit vector lengths, to get the actual PSP distances.
        dists_psp = (dists_single_dir[..., 0::2] + dists_single_dir[..., 1::2]) * len_dir_vectors[
                                                                                  :num_directions]
        # Select and count PSP events that are shorter than the given cutoff.
        psp_is_within_range = np.logical_and(has_psp, dists_psp <= psp_len_max)
        return np.count_nonzero(psp_is_within_range, axis=-1)



def routine_run_t2fpharm(
        receptor_filepath: Path,
        grid_center: Union[Tuple[float, float, float], Literal["auto"]],
        grid_size: Tuple[float, float, float] = (16, 16, 16),
        grid_spacing: float = 0.6,
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        param_filepath: Optional[Path] = None,
        output_path: Optional[Path] = None,
        vacancy_max_energy: float = +0.6,
        buriedness_psp_num_dirs: Literal[3, 7, 13] = 7,
        buriedness_psp_max_len: float = 10.0,
        buriedness_psp_min_count: int = 4,
        hbond_max_len: float = 3.0,
        probes_max_energy: Tuple[float, float, float, float] = (-0.6,-0.35,-0.4,-0.4),
        electrostat_pot_exclusion_range: Tuple[float, float] = (-1.0, 1.0),
        min_neighbor_dist_clustering: float = 1.21,
        min_common_neighbor_count_clustering: int = 6,
        min_points_per_cluster_count: int = 15,

):
    hba = autodock.AtomType.OA  # hydrogen-bond acceptor probe
    hbd = autodock.AtomType.HD  # hydrogen-bond donor probe
    aliph = autodock.AtomType.C  # aliphatic hydrophobic probe
    arom = autodock.AtomType.A  # aromatic hydrophobic probe
    probes = (hba, hbd, aliph, arom)

    t2fpharm = T2FPharm(
        receptor_filepath=receptor_filepath,
        ligand_types=probes,
        grid_center=grid_center,
        grid_size=grid_size,
        grid_spacing=grid_spacing,
        smooth=smooth,
        dielectric=dielectric,
        param_filepath=param_filepath,
        output_path=output_path,
    )
    # Vacancy of each grid point is a boolean grid
    is_vacant = t2fpharm.grid_probe(hba) < vacancy_max_energy
    # Buriedness of each vacant point as a boolean array
    is_buried = t2fpharm.count_psp_events(
        mask_vacancy=is_vacant,
        num_directions=buriedness_psp_num_dirs,
        psp_len_max=buriedness_psp_max_len,
    ) >= buriedness_psp_min_count
    # A boolean grid indicating points that are both vacant and buried
    is_vacant_n_buried = np.zeros_like(is_vacant)
    is_vacant_n_buried[is_vacant] = is_buried
    # A boolean grid indicating whether a specific probe at a specific grid point has high affinity
    is_high_affinity = t2fpharm.grid_probes < np.array(probes_max_energy)
    # A boolean grid indicating points that are both high affinity and buried
    is_high_affinity_n_buried = np.logical_and(
        is_high_affinity,
        is_vacant_n_buried[..., np.newaxis],
    )
    # Count of H-bond acceptor and donor atoms (from protein) in proximity of each grid point
    count_hba, count_hbd = t2fpharm.count_hbond_partners_in_target(
        max_len_hbond=hbond_max_len,
        mask_grid=is_vacant_n_buried
    )
    # Boolean grids describing the H-bond environment for each grid point
    has_hba_env = count_hba > 0
    has_hbd_env = count_hbd > 0
    has_hydrophilic_env = np.logical_or(has_hba_env, has_hbd_env)
    has_hydrophobic_env = np.logical_not(has_hydrophilic_env)

    is_hba, is_hbd, is_aliph, is_arom = (
        np.logical_and(is_high_affinity_n_buried[..., probes.index(probe)], condition)
        for probe, condition in zip(probes, [has_hbd_env, has_hba_env, has_hydrophobic_env, has_hydrophobic_env])
    )


    mask_pi = np.logical_and(
        is_vacant_n_buried,
        t2fpharm.grid_electrostat_pot < -elec_max_abs_pot,
        has_hydrophilic_env,
    )

    mask_ni = np.logical_and(
        is_vacant_n_buried,
        t2fpharm.grid_electrostat_pot > elec_max_abs_pot,
        has_hydrophilic_env
    )
