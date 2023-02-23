from typing import Optional, Union, Sequence, Literal
import opencadd as oc

import numpy as np
import jax.numpy as jnp


class LigSiteDetector:

    def __init__(
            self,
            receptor: oc.chem.system.ChemicalSystem,
            resolution_or_grid: Union[float, Sequence[float], oc.spacetime.grid.Grid],
            num_directions: Literal[3, 7, 13] = 13,
            max_radius: Optional[float] = None,
    ):
        self._receptor = receptor
        self._toxel_vol = self._receptor.conformation.toxelate(
            resolution_or_grid=resolution_or_grid,
            radius_points=1.7,
        )

        dimensions = {3: 1, 7: (1, 3), 13: None}
        dir_vectors = self._toxel_vol.spatial_direction_vectors(
            dimensions=dimensions[num_directions]
        )
        len_dir_vectors = np.linalg.norm(dir_vectors[:, 1:] * self._toxel_vol.grid.spacings, axis=-1)
        num_dir_vectors = dir_vectors.shape[0]
        # Calculate distance of each vacant grid point to the nearest occupied grid point in each
        # half direction, in units of corresponding distance vectors
        dists = self._toxel_vol.xeno_neighbor_distance(
            dir_vectors=dir_vectors,
            dir_multipliers=(max_radius // len_dir_vectors) if max_radius is not None else None
        ).astype(np.single)
        # set distances that are 0 (meaning no neighbor was found in that direction) to Nan.
        dists[dists == 0] = np.nan
        self._ps_dists_in_dir_units = dists
        self._ps_dists = dists * len_dir_vectors
        # Add distances to neighbors in positive half-directions , to distances to neighbors in
        # negative half-directions, in order to get the PSP length in units of direction vectors,
        # and then multiply by direction unit vector lengths, to get the actual PSP distances.
        self._psp_distances = dists[..., :num_dir_vectors // 2] + dists[..., num_dir_vectors // 2:]
        return

    @property
    def toxel_vol(self):
        return self._toxel_vol

    def routine_run(self):
        return oc.pocket.BindingPocket(volume=self.calculate_buriedness())

    def calculate_buriedness(
            self,
            psp_max_length: float = 20.0,
            psp_min_count: int = 2,
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
        # Count PSP events that are shorter than the given cutoff.
        grid_psp_counts = jnp.count_nonzero(self._psp_distances <= psp_max_length, axis=-1)
        buriedness = grid_psp_counts >= psp_min_count
        site = jnp.logical_and(jnp.logical_not(self.toxel_vol.toxels), buriedness)
        return site

    def psp_events(
            self,
            num_directions: Literal[3, 7, 13] = 13,
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

        return
