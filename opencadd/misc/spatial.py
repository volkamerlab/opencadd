"""
General functions and routines for spatial calculations, such as distances, volumes, etc.
"""


# Standard library
from typing import Any, Sequence, Union, Optional, Tuple
import itertools
# 3rd-party
import numpy as np
from scipy.spatial import distance_matrix


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


class Grid:
    def __init__(
            self,
            coords_origin: Sequence[float],
            shape: Sequence[int],
            spacing: float,
            name_properties: Tuple[Union[str, int]],
            dtype: np.dtype = np.single,
    ):
        if len(coords_origin) != len(shape):
            raise ValueError("`coords_origin` and `shape` must have the same dimension.")
        self._coords_origin = coords_origin
        self._shape = shape
        self._spacing = spacing
        if len(name_properties) > 1:
            self._name_properties = np.array(name_properties)
            tensor_shape = (*shape, len(name_properties))
        else:
            self._name_properties = np.array([name_properties])
            tensor_shape = shape
        self._tensor = np.zeros(
            shape=tensor_shape,
            dtype=dtype
        )

    @property
    def shape(self):
        return self._shape

    @property
    def spacing(self):
        return self._spacing

    @property
    def name_properties(self):
        return self._name_properties

    def get_properties(self, names: Tuple[Union[str, int]]):
        indices = np.where(np.isin(self._name_properties, names))[0]
        indices_squeezed = indices if len(indices) > 1 else indices[0]
        return self._tensor[..., indices_squeezed]

    def set_property(self, name, values):
        if np.isin(name, self._name_properties, invert=True):
            raise ValueError("Property was not found.")
        if len(self._name_properties) == 1:
            self._tensor[...] = values
        else:
            self._tensor[..., np.where(self._name_properties == name)[0][0]] = values
        return

    def coordinates_grid_points(
            self,
            mask: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        Calculate coordinates of the grid points, in the reference frame of the target structure.
        These are calculated from the same input parameters used in function `create_gpf` to create
        the configuration file for AutoGrid.

        Parameters
        ----------
        gridcenter : Sequence[float, float, float], Optional, default: "auto"
            Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
        npts : Sequence[int, int, int], Optional, default: (40, 40, 40)
            Number of grid points to add to the central grid point, along x-, y- and z-axes,
            respectively. Each value must be an even integer number; when added to the central grid
            point, there will be an odd number of points in each dimension. The number of x-, y and
            z-grid points need not be equal.
        spacing : float, Optional, default: 0.375
            The grid-point spacing, i.e. distance between two grid points in Ångstrom (Å).
            Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
            all three dimensions.

        Returns
        -------
        numpy.ndarray
            A 4-dimensional array containing the (x,y,z)-coordinates of each grid point
            in Ångstrom (Å), in the target structure's reference frame.

            The shape of the array is (nx, ny, nz, 3), where nx, ny, and nz are the number of grid
            points along the x-, y-, and z-axis, respectively. These are each equal to the
            corresponding input `npts` value, plus 1 (due to center point). The array is ordered in
            the same way that the grid points are ordered in space, and thus the individual grid
            points can be indexed using their actual coordinates (in unit vectors), assuming the
            origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
            grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.
        """
        if mask is None:
            mask = np.ones(shape=self._shape, dtype=np.bool_)
        return np.argwhere(mask) * self._spacing + self._coords_origin

    def calculate_mask_proximate(
            self,
            coords_targets: np.ndarray,
            distance: float,
            mask: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        dist_matrix = distance_matrix(
            self.coordinates_grid_points(mask=mask),
            coords_targets
        )
        return dist_matrix < distance

    def select_by_property(
            self,
            names: Tuple[Union[str, int]],

    ):

    def __getitem__(self, item):
        return self._tensor.__getitem__(item)