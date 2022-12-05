"""
T2F (Truly Target-Focused) Pharmacophore modeler.

This module contains the `T2Fpharm` class, used for pharmacophore modeling
from protein apo structures.
"""


# Standard library
from typing import Literal, Optional, Sequence, Tuple, Union
from pathlib import Path
# 3rd-party
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import nglview
# Self
from opencadd import protein
from opencadd.interaction import autogrid
from opencadd.consts import autodock
from opencadd.misc import spatial
from opencadd.visualization import nglview_api
from opencadd.typing import PathLike


class T2FPharm:
    """
    Truly Target Focused (T2F) pharmacophore modeler.

    Design a pharmacophore from one or several protein apo structures.

    Instantiate the modeler by entering one or several protein structures, defining their
    binding pocket's coordinates and size, and indicating a set of pharmacophore features
    to consider. The modeler automatically calculates interaction energy fields inside the
    binding pocket of each protein structure, for each pharmacophore feature, plus desolvation
    energy and electrostatic potential fields. All data are stored in a single higher-dimensional
    array, allowing for easy and fast data manipulation and visualization. With the methods
    available to the class, a variety of information on the binding pocket can be extracted and
    used to develop a pharmacophore model.



    calculate the interaction energies
    between each protein structure and each atomic probe (representing a pharmacophore feature),
    at regular grid points spanned over the defined binding pocket.

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
            receptor: protein.Protein,
            interaction_field: autogrid.IntraMolecularInteractionField,
    ):
        """
        Parameters
        ----------
        receptors : Sequence[pathlib.Path]
            Path to the PDBQT structure files of the macromolecule.
        probe_types : Sequence[opencadd.consts.autodock.AtomType]
            Types of ligand atoms, for which interaction energies must be calculated.
            For more information, see `opencadd.consts.autodock.AtomType`.
        receptor_pocket_center : tuple[float, float, float] | Literal["auto"], Optional, default: "auto"
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
        self._receptor: protein.Protein = receptor
        self._interaction_field = interaction_field

        self._vacancy: np.ndarray = None
        self._psp_distances: np.ndarray = None
        self._buriedness: np.ndarray = None
        self._distances_to_protein_atoms: np.ndarray = None
        self._hbonding_atoms_count: np.ndarray = None
        return

    @property
    def receptor(self) -> protein.Protein:
        return self._receptor

    @property
    def field(self):
        return self._interaction_field

    @property
    def distances_to_protein_atoms(self):
        return (
            self._distances_to_protein_atoms if self._distances_to_protein_atoms is not None
            else self.calculate_distances_to_protein_atoms()
        )

    @property
    def vacancy(self):
        return self._vacancy if self._vacancy is not None else self.calculate_vacancy()

    @property
    def psp_distances(self):
        return self._psp_distances if self._psp_distances is not None else self.calculate_psp_distances()

    @property
    def buriedness(self):
        return self._buriedness if self._buriedness is not None else self.calculate_buriedness()

    @property
    def hbonding_atoms_count(self):
        return self._hbonding_atoms_count if self._hbonding_atoms_count is not None else self.count_hbonding_atoms_in_radius()

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
        vacant = self.vacancy(energy_cutoff=vacancy_max_energy)
        # Buriedness of each vacant point as a boolean array
        buried = self.calculate_buriedness(
            vacancy=vacant,
            psp_distances=self.grid_psp_dist(
                grid_vacancy=vacant,
                num_directions=buriedness_psp_num_dirs,
                max_radius=buriedness_psp_max_len
            ),
            psp_max_length=buriedness_psp_max_len,
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
        count_hba, count_hbd = self.count_hbonding_atoms_in_radius(
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

    def calculate_vacancy(
            self,
            energy_cutoff: float = +0.6,
            mode: Optional[Literal["max", "min", "avg", "sum"]] = "min",
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
        # if ligand_types is None:
        #     ind = slice(None)
        # else:
        #     ind = np.argwhere(np.expand_dims(ligand_types, axis=-1) == self._probe_types)[:, 1]
        #     # Verify that all input ligand types are valid
        #     if len(ind) != len(ligand_types):
        #         raise ValueError(f"Some of input energies were not calculated.")
        # Reduce the given references using the given operation.
        energy_vals = red_fun[mode](self._interaction_field.van_der_waals, axis=-1)
        # Apply cutoff and return
        self._vacancy = energy_vals < energy_cutoff
        return self._vacancy

    def calculate_psp_distances(
            self,
            vacancy: Optional[np.ndarray] = None,
            num_directions: Literal[3, 7, 13] = 7,
            max_radius: Optional[float] = None,
    ) -> np.ndarray:
        """
        Calculate the length of protein–solvent–protein (PSP) events for vacant grid points,
        and the length of solvent–protein–solvent events for occupied grid points,
        in each direction.

        Parameters
        ----------
        vacancy : numpy.ndarray
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
        if vacancy is None:
            vacancy = self.vacancy
        dimensions = {3: 1, 7: (1, 3), 13: None}
        dir_vectors = self._interaction_field.spatial_direction_vectors(
            dimensions=dimensions[num_directions]
        )
        len_dir_vectors = np.linalg.norm(dir_vectors, axis=-1)
        num_dir_vectors = dir_vectors.shape[0]
        # Calculate distance of each vacant grid point to the nearest occupied grid point in each
        # half direction, in units of corresponding distance vectors
        dists = spatial.xeno_neighbor_distance(
            bool_array=vacancy,
            dir_vectors=dir_vectors,
            dir_multipliers=(max_radius // len_dir_vectors) if max_radius is not None else None
        ).astype(np.single)
        # set distances that are 0 (meaning no neighbor was found in that direction) to Nan.
        dists[dists == 0] = np.nan
        # Add distances to neighbors in positive half-directions , to distances to neighbors in
        # negative half-directions, in order to get the PSP length in units of direction vectors,
        # and then multiply by direction unit vector lengths, to get the actual PSP distances.
        psp_grid_dists = dists[..., :num_dir_vectors//2] + dists[..., num_dir_vectors//2:]
        self._psp_distances = psp_grid_dists * len_dir_vectors[:num_dir_vectors//2]
        return self._psp_distances

    def calculate_buriedness(
            self,
            vacancy: Optional[np.ndarray] = None,
            psp_distances: Optional[np.ndarray] = None,
            psp_max_length: float = 10.0,
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
        vacancy : numpy.ndarray

        psp_distances : numpy.ndarray
        psp_max_length : float, Optional, default: 10.0
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
        if vacancy is None:
            vacancy = self.vacancy
        if psp_distances is None:
            psp_distances = self.psp_distances
        # Count PSP events that are shorter than the given cutoff.
        grid_psp_counts = np.count_nonzero(psp_distances <= psp_max_length, axis=-1)
        buriedness = grid_psp_counts >= psp_min_count
        buriedness[np.logical_not(vacancy)] = False
        return buriedness

    def calculate_distances_to_protein_atoms(self):
        self._distances_to_protein_atoms = self._interaction_field.grid.distance(
            coordinates=self._receptor.trajectory
        )
        return self._distances_to_protein_atoms

    def count_hbonding_atoms_in_radius(
            self,
            max_len_hbond: float = 3.0,
    ) -> np.ndarray:
        proximates = self.distances_to_protein_atoms < max_len_hbond
        self._hbonding_atoms_count = np.zeros(
            shape=(*self.distances_to_protein_atoms.shape[:-1], 2),
            dtype=np.byte
        )
        for idx, hbond_role in enumerate(["hbond_acc", "hbond_don"]):
            self._hbonding_atoms_count[..., idx] = np.count_nonzero(
                np.logical_and(self._receptor.atom_data[hbond_role], proximates),
                axis=-1
            )
        return self._hbonding_atoms_count

    def visualize(
            self,
            grid_mask = None,
            weights1 = None,
            weights2 = None,
            view = None,
            color_map: str = "bwr",
            opacity: float = 0.8,
    ):
        if grid_mask is None:
            grid_mask = np.ones(shape=self._interaction_field.grid.shape, dtype=np.bool_)
        if view is None:
            view = self.receptor.view
        normalizer = mpl.colors.Normalize()
        mapper = plt.cm.ScalarMappable(
            norm=normalizer,
            cmap=mpl.colormaps[color_map]
        )
        if weights1 is None:
            weights1 = np.repeat(0.5, np.count_nonzero(grid_mask) * 3)
        elif isinstance(weights1, Sequence):
            weights1 = np.tile(weights1, np.count_nonzero(grid_mask) * 3)
        else:
            weights1 = (mapper.to_rgba(weights1.flatten())[..., :3]).flatten()
        if weights2 is None:
            weights2 = (
                    np.ones(np.count_nonzero(grid_mask)) * self._interaction_field.grid.spacing / 4
            )
        elif isinstance(weights2, (int, float)):
            weights2 = np.ones(np.count_nonzero(grid_mask)) + weights2
        else:
            weights2 = normalizer(weights2).flatten() * (self._interaction_field.grid.spacing - 0.05) + 0.05
        nglview_api.add_spheres(
            view=view,
            coords=list(self._interaction_field.grid.coordinates[grid_mask].flatten()),
            colors=list(weights1),
            radii=list(weights2),
            opacity=opacity
        )
        return view

# RECEPTOR_FILEPATHS = [Path("/Users/home/Downloads/3w32.pdbqt")]
# GRID_CENTER = (15.91, 32.33, 11.03)
# GRID_SIZE = (20, 20, 20)
# GRID_SPACING = 0.6
# HBD = autodock.AtomType.HD
# HBA = autodock.AtomType.OA
# AL = autodock.AtomType.C
# AR = autodock.AtomType.A
# PROBES = (HBD, HBA, AL, AR)
# VACANCY_MAX_ENERGY = 0.6
# PSP_NUM_DIRS = 7
# PSP_LEN_MAX = 10.0
# PSP_COUNT_MIN = 4
# HBOND_LEN_MAX = 3.0
# t2f = T2FPharm(
#     receptor_filepaths=RECEPTOR_FILEPATHS,
#     ligand_types=PROBES,
#     grid_center=GRID_CENTER,
#     grid_size=GRID_SIZE,
#     grid_spacing=GRID_SPACING
# )
# t2f()



        #
        #
        # if isinstance(receptors, (PathLike, Sequence)):
        #     if isinstance(receptors, PathLike):
        #         path = Path(receptors)
        #         if path.is_file():
        #             if path.suffix.lower() == ".pdbqt":
        #                 pdbqt_filepaths = [path]
        #             else:
        #                 raise ValueError("Receptor's input file's extension should be PDBQT.")
        #         elif path.is_dir():
        #             pdbqt_filepaths = list(path.glob("*.pdbqt"))
        #         else:
        #             raise ValueError(f"No such file or directory: {path}")
        #     else:
        #         pdbqt_filepaths = [
        #             Path(receptor) for receptor in receptors
        #             if Path(receptor).is_file() and Path(receptor).suffix.lower() == ".pdbqt"
        #         ]
        #     if len(pdbqt_filepaths) == 0:
        #         raise ValueError(f"No PDBQT file found.")
        #     if len(pdbqt_filepaths) == 1:
        #         self._receptor = protein.Protein.from_file(path=pdbqt_filepaths[0])
        #         self._receptor_is_static = True
        #         self._receptor_is_trajectory = False
        #     else:
        #         self._receptor_is_static = False
        #         try:
        #             self._receptor = protein.ProteinTrajectory.from_files(paths=pdbqt_filepaths)
        #             self._receptor_is_trajectory: bool = True
        #         except protein.TopologyError:
        #             self._receptor = protein.ProteinFamily.from_files(paths=pdbqt_filepaths)
        #             self._receptor_is_trajectory: bool = False
        # elif isinstance(
        #         receptors,
        #         (protein.Protein, protein.ProteinTrajectory, protein.ProteinFamily)
        # ):
        #     self._receptor = receptors
        #     pdbqt_filepaths = self._receptor.to_file_pdbqt()
        #     if isinstance(receptors, protein.ProteinFamily):
        #         self._receptor_is_static = False
        #         self._receptor_is_trajectory = False
        #     elif isinstance(receptors, protein.ProteinTrajectory):
        #         self._receptor_is_static = False
        #         self._receptor_is_trajectory = True
        #     else:
        #         self._receptor_is_static = True
        #         self._receptor_is_trajectory = False
        #         pdbqt_filepaths = [pdbqt_filepaths]
        # else:
        #     raise ValueError("Type of `receptor` not recognized.")
        #
        #
        # self._output_path = pdbqt_filepaths[0].parent if output_path is None else output_path
        # self._output_path.mkdir(parents=True, exist_ok=True)
        # self._probe_types = np.array(probe_types, dtype=object)
        #
        #
        #
        #
        # self._df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptors[0])["ATOM"]
        #
        # for t, receptor_filepath in enumerate(receptors):
        #     df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptor_filepath)["ATOM"]
        #     self._protein_coordinates[t, ...] = df_pdbqt[["x", "y", "z"]].to_numpy()
        #
        # self._grid_dists_protein = np.moveaxis(
        #     distance_matrix(
        #         self._grid_coords.reshape(-1, 3),
        #         self._protein_coordinates.reshape(-1, 3)
        #     ).reshape(*self._grid_shape, self._temporal_len, len(self._df_pdbqt.index)),
        #     source=3,
        #     destination=0,
        # )
        #
        # # Get all 26 half directions, in the order: orthogonal, 3d-diagonal, 2d-diagonal
        # self._direction_vects = np.zeros(shape=(26, 4), dtype=np.byte)
        # self._direction_vects[:, 1:] = np.concatenate(
        #     [spatial.GRID_DIRS_ORTHO, spatial.GRID_DIRS_DIAG_3D, spatial.GRID_DIRS_DIAG_2D]
        # )
        # # Calculate the length of each direction unit vector
        # self._direction_vects_len = np.linalg.norm(self._direction_vects, axis=-1) * self._grid_spacing