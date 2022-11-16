"""
General functions and routines for spatial calculations, such as distances, volumes, etc.
"""


# Standard library
from typing import Any, Sequence
# 3rd-party
import numpy as np


# Direction vectors in a 3x3x3 grid:
GRID_DIR_ORTHO_POS = np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.byte)
GRID_DIR_DIAG_2D_POS = np.array(
    [[1,1,0], [1,0,1], [0,1,1], [-1,1,0], [-1,0,1], [0,-1,1]],
    dtype=np.byte
)
GRID_DIR_DIAG_3D_POS = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=np.byte)
GRID_DIRS_ORTHO_NEG = -GRID_DIR_ORTHO_POS
GRID_DIRS_DIAG_2D_NEG = -GRID_DIR_DIAG_2D_POS
GRID_DIRS_DIAG_3D_NEG = -GRID_DIR_DIAG_3D_POS
GRID_DIRS_ORTHO = np.insert(GRID_DIRS_ORTHO_NEG, np.arange(3), GRID_DIR_ORTHO_POS, axis=0)
GRID_DIRS_DIAG_2D = np.insert(GRID_DIRS_DIAG_2D_NEG, np.arange(6), GRID_DIR_DIAG_2D_POS, axis=0)
GRID_DIRS_DIAG_3D = np.insert(GRID_DIRS_DIAG_3D_NEG, np.arange(4), GRID_DIR_DIAG_3D_POS, axis=0)
GRID_DIRS_DIAG = np.concatenate([GRID_DIRS_DIAG_2D, GRID_DIRS_DIAG_3D])
GRID_DIRS = np.concatenate([GRID_DIRS_ORTHO, GRID_DIRS_DIAG])


def grid_distance(
        grid: np.ndarray,
        start_indices: np.ndarray,
        directions: np.ndarray,
        target_value: Any,
        max_dist_in_dir: Sequence[int],
) -> np.ndarray:
    """
    On an n-dimensional grid of values, for all given starting grid points, find the grid
    distance (i.e. number of intermediate grid points plus 1) to the nearest point with a given
    target value, in each of the given directions, optionally up to a maximum grid distance.

    Parameters
    ----------
    grid : numpy.ndarray
        An n-dimenstional array representing the values of the grid points.
    start_indices : numpy.ndarray[dtype=int]
        A 2-dimensional array of shape (i, n), consisting of coordinates (set of indices) of 'i'
        gird points, for which the distance to nearest neighbors are to be found.
    directions : numpy.ndarray[dtype=int]
        A 2-dimensional array of shape (k, n), consisting of 'k' n-dimensional direction
        vectors, along which the nearest target grid-points are to be found. The direction vectors
        indicate the number of grid points to travers in each dimension to find the first
        neighboring point in that direction, and must be (positive or negative) integers.
        For example in 3d space, [1,0,0] is along the positive x-axis and [0,-1, 0] is along the
        negative y-axis, whereas [1,1,0] goes to the next diagonal point along positive xy-plane.
    target_value : Any
        The value of the grid points, to which the distances are to be calculated.
    max_dist_in_dir : Sequence[int]
        A sequence of length 'k', indicating the maximum number of neighboring points to search, in
        each direction in `directions`.

    Returns
    -------
    dists : numpy.ndarray[dtype=numpy.ushort]
        A 2-dimensional array of shape (i, k), where each element at position (m, n) represents
        the minimum number of times to move from the 'm'th point in `start_indices` along the
        'n'th direction vector in `directions`, in order to land on a point with value
        `target_value`. A value of 0 means that no point with the given target value is found in
        that direction, because either the end of the grid or `max_dist_in_dir` is reached.
    """
    # TODO: vectorize this
    dists = np.zeros((start_indices.shape[0], directions.shape[0]), dtype=np.uintc)
    for idx_point, point in enumerate(start_indices):
        for idx_dir, direct in enumerate(directions):
            len_max = np.min(
                np.concatenate([(np.array(grid.shape) - 1 - point)[direct==1], point[direct==-1]])
            )
            for i in range(1, int(min(len_max, max_dist_in_dir[idx_dir]))+1):
                if grid[tuple(point + i * direct)] == target_value:
                    dists[idx_point, idx_dir] = i
                    break
    return dists
