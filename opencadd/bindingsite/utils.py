"""
This module contains general utility functions used in this subpackage.
"""


# Standard library
from typing import Literal
# 3rd-party
import numpy as np
# Self
from opencadd.misc.spatial import grid_distance, GRID_DIRS_ORTHO, GRID_DIRS_DIAG_2D, GRID_DIRS_DIAG_3D


def count_psp_events_from_vacancy_tensor(
            grid_vacancy: np.ndarray,
            spacing_grid_points: float,
            num_directions: Literal[3, 7, 13] = 7,
            psp_len_max: float = 10.0,
    ) -> np.ndarray:
    """
    Count the number of protein-solvent-protein (PSP) events for each point on a 3-dimensional
    grid, given a 3-dimensional boolean array indicating the vacancy of each grid point.

    Parameters
    ----------
    grid_vacancy : numpy.ndarray
        A 3-dimensional boolean array representing a 3-dimensional grid, indicating
        whether each grid point is vacant (True), or occupied (False).
    spacing_grid_points : float
        The grid-point spacing, i.e. the distance between two adjacent grid points.
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
    dirs = np.concatenate([GRID_DIRS_ORTHO, GRID_DIRS_DIAG_3D, GRID_DIRS_DIAG_2D])
    # Calculate the length of each direction unit vector
    len_dir_vectors = np.repeat(spacing_grid_points, 26) * np.repeat(
        np.array([1, np.sqrt(3), np.sqrt(2)]),
        repeats=[6, 8, 12],
    )
    # Calculate the maximum allowed number of times to travel in each half-direction
    max_mult_factors = psp_len_max // len_dir_vectors
    # Calculate distance of each vacant grid point to the nearest occupied grid point in each
    # half direction, in units of corresponding distance vectors
    dists_single_dir = grid_distance(
        grid=grid_vacancy,
        start_indices=np.argwhere(grid_vacancy),
        directions=dirs[:num_directions*2],
        target_value=0,
        max_dist_in_dir=max_mult_factors[:num_directions*2]
    )
    # For each grid-point, filter directions where one half-direction is zero (meaning no
    # neighbor was found in that direction).
    has_psp = np.logical_and(dists_single_dir[..., 0::2] != 0, dists_single_dir[..., 1::2] != 0)
    # Add distances to neighbors in positive half-directions , to distances to neighbors in
    # negative half-directions, in order to get the PSP length in units of direction vectors,
    # and then multiply by direction unit vector lengths, to get the actual PSP distances.
    dists_psp = (dists_single_dir[...,0::2] + dists_single_dir[...,1::2]) * len_dir_vectors[:num_directions]
    # Select and count PSP events that are shorter than the given cutoff.
    psp_is_within_range = np.logical_and(has_psp, dists_psp <= psp_len_max)
    return np.count_nonzero(psp_is_within_range, axis=-1)


