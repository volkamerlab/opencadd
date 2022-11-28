"""
T2F-Pharm (Truly Target-Focused Pharmacophore) model.
"""


# Standard library
from typing import Sequence, Union, Optional, Literal, Tuple, Dict
from pathlib import Path
# 3rd-party
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import nglview
from scipy.spatial import distance_matrix
# Self
from opencadd.docking import autogrid
from opencadd.io import pdbqt
from opencadd.consts import autodock
from opencadd.misc import spatial
from opencadd.visualization import nglview_api


class T2FPharm:
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
            receptor_filepaths: Sequence[Path],
            ligand_types: Sequence[autodock.AtomType],
            grid_center: Union[Tuple[float, float, float], Literal["auto"]],
            grid_size: Tuple[float, float, float] = (16, 16, 16),
            grid_spacing: float = 0.6,
            smooth: float = 0.5,
            dielectric: float = -0.1465,
            param_filepath: Optional[Path] = None,
            output_path: Optional[Path] = None,
            data_type: np.dtype = np.single,
    ):
        """

        Parameters
        ----------
        receptor_filepaths : Sequence[pathlib.Path]
            Path to the PDBQT structure files of the macromolecule.
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
        self._receptor_filepaths = receptor_filepaths
        self._temporal_len = len(receptor_filepaths)
        self._output_path = receptor_filepaths[0].parent if output_path is None else output_path
        self._output_path.mkdir(parents=True, exist_ok=True)
        self._ligand_types = np.array(ligand_types, dtype=object)
        self._grid_center = np.array(grid_center)
        self._grid_size = grid_size
        self._grid_spacing = grid_spacing
        npts = np.array(autogrid.calculate_npts(grid_size=grid_size, grid_spacing=grid_spacing))
        self._grid_shape = npts + 1
        self._grid_origin = self._grid_center - self._grid_spacing * npts / 2
        self._grid = np.empty(
            shape=(self._temporal_len, *self._grid_shape, self._ligand_types.size+2),
            dtype=data_type
        )
        self._df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptor_filepaths[0])["ATOM"]
        self._protein_coordinates = np.empty(
            shape=(self._temporal_len, len(self._df_pdbqt.index), 3)
        )
        for t, receptor_filepath in enumerate(receptor_filepaths):
            df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptor_filepath)["ATOM"]
            self._protein_coordinates[t, ...] = df_pdbqt[["x", "y", "z"]].to_numpy()
            energies, filepaths_mapfiles = autogrid.routine_run(
                receptor_filepath=receptor_filepath,
                ligand_types=ligand_types,
                grid_center=grid_center,
                grid_npts=tuple(npts),
                grid_spacing=grid_spacing,
                smooth=smooth,
                dielectric=dielectric,
                param_filepath=param_filepath,
                output_path=output_path,
            )
            for idx_energy, energy in enumerate(energies):
                self._grid[t, ..., idx_energy] = energy.reshape(tuple(self._grid_shape), order="F")

        self._grid_coords = np.array(
            list(np.ndindex(*self._grid_shape))
        ).reshape(*self._grid_shape, 3) * self._grid_spacing + self._grid_origin

        self._grid_dists_protein = np.moveaxis(
            distance_matrix(
                self._grid_coords.reshape(-1, 3),
                self._protein_coordinates.reshape(-1, 3)
            ).reshape(*self._grid_shape, self._temporal_len, len(self._df_pdbqt.index)),
            source=3,
            destination=0,
        )

        # Get all 26 half directions, in the order: orthogonal, 3d-diagonal, 2d-diagonal
        self._direction_vects = np.zeros(shape=(26, 4), dtype=np.byte)
        self._direction_vects[:, 1:] = np.concatenate(
            [spatial.GRID_DIRS_ORTHO, spatial.GRID_DIRS_DIAG_3D, spatial.GRID_DIRS_DIAG_2D]
        )
        # Calculate the length of each direction unit vector
        self._direction_vects_len = np.linalg.norm(self._direction_vects, axis=-1) * self._grid_spacing
        return

    def __str__(self):
        str = f"""
        Grid shape (n_, nx, ny, nz, n)
        """

    def __call__(
            self,
            vacancy_max_energy: float = +0.6,
            buriedness_psp_num_dirs: Literal[3, 7, 13] = 7,
            buriedness_psp_max_len: float = 10.0,
            buriedness_psp_min_count: int = 4,
            hbond_max_len: float = 3.0,
            probes_max_energy: Tuple[float, float, float, float] = (-0.6, -0.35, -0.4, -0.4),
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

        # Vacancy of each grid point as a boolean grid
        vacant = self.grid_vacancy(energy_cutoff=vacancy_max_energy)
        # Buriedness of each vacant point as a boolean array
        buried = self.grid_buriedness(
            grid_vacancy=vacant,
            grid_psp_dist=self.grid_psp_dist(
                grid_vacancy=vacant,
                num_directions=buriedness_psp_num_dirs,
                max_radius=buriedness_psp_max_len
            ),
            psp_max_len=buriedness_psp_max_len,
            psp_min_count=buriedness_psp_min_count
        )
        # A boolean grid indicating whether a specific probe at a specific grid point has high affinity
        is_high_affinity = self.grid_probes < np.array(probes_max_energy)
        # A boolean grid indicating points that are both high affinity and buried
        high_affinity_n_buried = np.logical_and(
            is_high_affinity,
            buried[..., np.newaxis],
        )
        # Count of H-bond acceptor and donor atoms (from protein) in proximity of each grid point
        count_hba, count_hbd = self.grid_hbonding_atoms_count(
            max_len_hbond=hbond_max_len,
        )
        # Boolean grids describing the H-bond environment for each grid point
        has_hba_env = count_hba > 0
        has_hbd_env = count_hbd > 0
        has_hydrophilic_env = np.logical_or(has_hba_env, has_hbd_env)
        has_hydrophobic_env = np.logical_not(has_hydrophilic_env)

        is_hba, is_hbd, is_aliph, is_arom = (
            np.logical_and(high_affinity_n_buried[..., probes.index(probe)], condition)
            for probe, condition in
        zip(probes, [has_hbd_env, has_hba_env, has_hydrophobic_env, has_hydrophobic_env])
        )

        is_buried_in_hydrophilic_env = np.logical_and(buried, has_hydrophilic_env)
        is_pi = np.logical_and(
            is_buried_in_hydrophilic_env,
            self.grid_electrostatic < electrostat_pot_exclusion_range[0],
        )
        is_ni = np.logical_and(
            is_buried_in_hydrophilic_env,
            self.grid_electrostatic > electrostat_pot_exclusion_range[1],
        )
        return

    @property
    def grid(self) -> np.ndarray:
        """
        5-dimensional array containing data of all grid points,
        calculated for each protein structure.

        Returns
        -------
        numpy.ndarray
            A 5-dimensional array of shape (n_t, n_x, n_y, n_z, n_p), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
            n_p: number of calculated properties for each point, i.e. number of input ligand types
                 plus 2 (for electrostatic potential and desolvation energy).
        """
        return self._grid

    @property
    def grid_probes(self) -> np.ndarray:
        """
        Sub-array of `T2FPharm.grid`, containing only the ligand interaction energies.

        Returns
        -------
        numpy.ndarray
            A 5-dimensional array of shape (n_t, n_x, n_y, n_z, n_l), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
            n_l: number of input ligand types.
        """
        return self._grid[..., :-2]

    @property
    def grid_electrostatic(self) -> np.ndarray:
        """
        Sub-array of `T2FPharm.grid`, containing only the electrostatic potential.

        Returns
        -------
        numpy.ndarray
            A 4-dimensional array of shape (n_t, n_x, n_y, n_z), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
        """
        return self._grid[..., -2]

    @property
    def grid_desolvation(self) -> np.ndarray:
        """
        Sub-array of `T2FPharm.grid`, containing only the desolvation energy.

        Returns
        -------
        numpy.ndarray
            A 4-dimensional array of shape (n_t, n_x, n_y, n_z), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
        """
        return self._grid[..., -1]

    @property
    def grid_coordinates(self):
        return self._grid_coords

    @property
    def receptor_data(self) -> pd.DataFrame:
        return self._df_pdbqt

    def grid_probe(self, ligand_type: autodock.AtomType) -> np.ndarray:
        """
        Sub-array of `T2FPharm.grid`, containing only
        the interaction energy of a given ligand type.

        Parameters
        ----------
        ligand_type : opencadd.consts.autodock.AtomType

        Returns
        -------
        numpy.ndarray
            A 4-dimensional array of shape (n_t, n_x, n_y, n_z), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
        """
        idx = np.argwhere(self._ligand_types == ligand_type)
        if idx.size == 0:
            raise IndexError("No energy has been calculated for the given ligand type.")
        return self._grid[..., idx[0, 0]]

    def grid_vacancy(
            self,
            energy_cutoff: float = +0.6,
            mode: Optional[Literal["max", "min", "avg", "sum"]] = "min",
            ligand_types: Optional[Sequence[autodock.AtomType]] = None,
    ) -> np.ndarray:
        """
        Calculate whether each grid point is vacant, or occupied by a target atom.

        Parameters
        ----------
        energy_cutoff : float, Optional, default: +0.6
            Cutoff value for energy; grid points with energies lower than cutoff are considered
            vacant.
        mode: Literal["max", "min", "avg", "sum"], Optional, default: "min"
            If the energy of more than one ligand type is to be compared, this parameter defines
            how those different energy values must be processed, before comparing with the cutoff.
        ligand_types : Sequence[opencadd.consts.autodock.AtomType], Optional, default: None
            A subset of ligand types that were used to initialize the object, whose energy values
            are to be taken as reference for calculating the vacancy of each grid point. If not
            set to None, then all ligand interaction energies are considered.

        Returns
        -------
        vacancy : numpy.ndarray[dtype=numpy.bool_, shape=T2FPharm.grid.shape[:-1]]
            A 4-dimensional boolean array matching the first four dimensions of `T2FPharm.grid`,
            indicating whether each grid point is vacant (True), or occupied (False).
            Vacant grid points can easily be indexed by `T2FPharm.grid[vacancy]`.
        """
        # The reducing operations corresponding to each `mode`:
        red_fun = {"max": np.max, "min": np.min, "avg": np.mean, "sum": np.sum}
        # Get index of input ligand types
        if ligand_types is None:
            ind = slice(None)
        else:
            ind = np.argwhere(np.expand_dims(ligand_types, axis=-1) == self._ligand_types)[:, 1]
            # Verify that all input ligand types are valid
            if len(ind) != len(ligand_types):
                raise ValueError(f"Some of input energies were not calculated.")
        # Reduce the given references using the given operation.
        energy_vals = red_fun[mode](self.grid_probes[..., ind], axis=-1)
        # Apply cutoff and return
        return energy_vals < energy_cutoff

    def grid_psp_dist(
            self,
            grid_vacancy: np.ndarray,
            num_directions: Literal[3, 7, 13] = 7,
            max_radius: Optional[float] = None,
    ) -> np.ndarray:
        """
        Calculate the length of protein–solvent–protein (PSP) events for vacant grid points,
        and the length of solvent–protein–solvent events for occupied grid points,
        in each direction.

        Parameters
        ----------
        grid_vacancy : numpy.ndarray
            A 4-dimensional boolean array, indicating whether each grid point is vacant (True),
            or occupied (False), i.e. the output of `T2FPharm.grid_vacancy`.
        num_directions : Literal[3, 7, 13]
            Number of directions to look for events for each grid point. Directions are
            all symmetric around each point, and each direction entails both positive and
            negative directions, along an axis. Directions are defined in a 3x3x3 unit cell,
            i.e. for each point, there are 26 neighbors, and thus max. 13 directions.
            Options are:
            3: x-, y- and z-directions
            7: x, y, z, and four 3d-diagonals (i.e. with absolute direction unit vector [1, 1, 1])
            13: x, y, z, four 3d-diagonals, and six 2d-diagonals (e.g. [1, 1, 0])
        max_radius : float, Optional, default: None
            Maximum radius of search, in the same unit as `spacing_grid_points`.

        Returns
        -------
        psp_dist : numpy.ndarray[dtype=numpy.single, shape=(*grid_vacancy.shape, num_directions)]
            A 4-dimensional array matching the first four dimensions of `T2FPharm.grid`,
            indicating for each grid point its psp distances in each given direction.
        """
        if num_directions not in [3, 7, 13]:
            raise ValueError("Argument `num_directions` must be either 3, 7, or 13.")
        # Calculate distance of each vacant grid point to the nearest occupied grid point in each
        # half direction, in units of corresponding distance vectors
        dists = spatial.dist_vectorized(
            grid=grid_vacancy,
            directions=self._direction_vects[:num_directions * 2],
            directions_mult=(
                    max_radius // self._direction_vects_len[:num_directions * 2]
            ) if max_radius is not None else None
        ).astype(np.single)
        # set distances that are 0 (meaning no neighbor was found in that direction) to Nan.
        dists[dists == 0] = np.nan
        # Add distances to neighbors in positive half-directions , to distances to neighbors in
        # negative half-directions, in order to get the PSP length in units of direction vectors,
        # and then multiply by direction unit vector lengths, to get the actual PSP distances.
        psp_grid_dists = dists[..., 0::2] + dists[..., 1::2]
        return psp_grid_dists * self._direction_vects_len[:num_directions * 2:2]

    @staticmethod
    def grid_buriedness(
            grid_vacancy: np.ndarray,
            grid_psp_dist: np.ndarray,
            psp_max_len: float = 10.0,
            psp_min_count: int = 4,
    ) -> np.ndarray:
        """
        Calculate whether each grid point is buried inside the target structure or not, based on
        counting the number of protein-solvent-protein (PSP) events for each point, and applying
        a cutoff.

        Notice that before using this method, a vacancy grid must be calculated using the method
        `calculate_vacancy_from_energy`.

        Parameters
        ----------
        grid_vacancy : numpy.ndarray

        grid_psp_dist : numpy.ndarray
        psp_max_len : float, Optional, default: 10.0
            Maximum acceptable distance for a PSP event, in Ångstrom (Å).
        psp_min_count : int, Optional, default: 4
            Minimum required number of PSP events for a grid point, in order to count as buried.

        Returns
        -------
        buriedness : numpy.ndarray[dtype=numpy.bool_, shape=grid_vacancy.shape]
            A 4-dimensional array matching the dimensions of the grid,
            indicating whether each grid point is buried (True), or exposed (False).
            Thus, the buried grid points can easily be indexed using boolean indexing with this
            array, i.e.: `grid[buriedness]`.
        """
        # Count PSP events that are shorter than the given cutoff.
        grid_psp_counts = np.count_nonzero(grid_psp_dist <= psp_max_len, axis=-1)
        buriedness = grid_psp_counts >= psp_min_count
        buriedness[np.logical_not(grid_vacancy)] = False
        return buriedness

    def grid_hbonding_atoms_count(
            self,
            max_len_hbond: float = 3.0,
    ) -> np.ndarray:
        proximates = self._grid_dists_protein < max_len_hbond
        counts = np.zeros(shape=(*self._grid.shape[:-1], 2), dtype=np.byte)
        for idx, hbond_role in enumerate(["hbond_acc", "hbond_don"]):
            counts[..., idx] = np.count_nonzero(
                np.logical_and(self._df_pdbqt[hbond_role], proximates),
                axis=-1
            )
        return counts


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
