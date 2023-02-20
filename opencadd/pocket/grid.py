from typing import Optional
from opencadd import Protein

import numpy as np


def embed_on_grid(

)


class PSPEvent:

    def __init__(
            self,
            vacancy_field,
            num_cubic_expansions: int = 1,
    ):
        pass


def calculate_protein_solvent_distances(
        vacancy_field: np.ndarray,
        direction_unit_vectors: np.ndarray,
        grid_spacing: float = 1,
        max_distance: Optional[float] = None,
) -> np.ndarray:
    pass

def psp_events(
        vacancy: BinaryToxelField,
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
    dir_vectors = self._fields.spatial_direction_vectors(
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



class GridBindingSiteDetector:

    def __init__(self, receptor: Protein):
        self._receptor: Protein = receptor
        return

