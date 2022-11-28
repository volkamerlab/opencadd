"""
General functions and routines for spatial calculations, such as distances, volumes, etc.
"""


# Standard library
from typing import Any, Sequence, Union, Optional, Tuple
import itertools
# 3rd-party
import numpy as np
from scipy.spatial import distance_matrix


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


def dist_vectorized(
        grid: np.ndarray,
        directions: np.ndarray,
        directions_mult: Optional[Union[Sequence[int], int]] = None,
):
    def slicer(vec):
        slice_start, slice_end, slice_exclude = [], [], []
        for val in vec:
            if val > 0:
                slice_start.append(slice(None, -val))
                slice_end.append(slice(val, None))
                slice_exclude.append(slice(-val, None))
            elif val < 0:
                slice_start.append(slice(-val, None))
                slice_end.append(slice(None, val))
                slice_exclude.append(slice(None, -val))
            else:
                for lis in [slice_start, slice_end, slice_exclude]:
                    lis.append(slice(None, None))
        return tuple(slice_start), tuple(slice_end), tuple(slice_exclude)

    dists = np.zeros(shape=(*grid.shape, directions.shape[0]), dtype=np.uintc)
    with np.errstate(divide='ignore', invalid='ignore'):
        max_mult_axis = (np.array(grid.shape) - 1) / np.abs(directions)
        max_mult_axis[np.isnan(max_mult_axis)] = np.inf
        max_mult_dir = np.min(max_mult_axis, axis=-1)
    if directions_mult is None:
        directions_mult = np.ones(directions.shape[0]) * np.max(grid.shape)
    elif isinstance(directions_mult, int):
        directions_mult = np.ones(directions.shape[0]) * directions_mult
    max_mult = np.min((max_mult_dir, directions_mult), axis=0).astype(int) + 1
    for idx_dir, direction in enumerate(directions):
        curr_mask = np.ones_like(grid)
        for mult in range(1, max_mult[idx_dir]):
            start_slice, end_slice, excl_slice = slicer(mult*direction)
            reached_xeno = np.logical_xor(grid[start_slice], grid[end_slice])
            dists[(*start_slice, idx_dir)][curr_mask[start_slice]] = reached_xeno[curr_mask[
                start_slice]] * mult
            curr_mask[start_slice][reached_xeno] = 0
    return dists



class Grid:
    """
    An n-dimensional grid of evenly spaced points.
    """
    # Direction vectors in a 3x3x3 grid:
    _DIR_ORTHO_POS = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.byte)
    _DIR_DIAG_2D_POS = np.array(
        [[1, 1, 0], [1, 0, 1], [0, 1, 1], [-1, 1, 0], [-1, 0, 1], [0, -1, 1]],
        dtype=np.byte
    )
    _DIR_DIAG_3D_POS = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]], dtype=np.byte)
    _DIRS_ORTHO_NEG = -_DIR_ORTHO_POS
    _DIRS_DIAG_2D_NEG = -_DIR_DIAG_2D_POS
    _DIRS_DIAG_3D_NEG = -_DIR_DIAG_3D_POS
    _DIRS_ORTHO = np.insert(_DIRS_ORTHO_NEG, np.arange(3), _DIR_ORTHO_POS, axis=0)
    _DIRS_DIAG_2D = np.insert(_DIRS_DIAG_2D_NEG, np.arange(6), _DIR_DIAG_2D_POS, axis=0)
    _DIRS_DIAG_3D = np.insert(_DIRS_DIAG_3D_NEG, np.arange(4), _DIR_DIAG_3D_POS, axis=0)
    _DIRS_DIAG = np.concatenate([_DIRS_DIAG_2D, _DIRS_DIAG_3D])
    _DIRS = np.concatenate([_DIRS_ORTHO, _DIRS_DIAG])

    def __init__(
            self,
            shape: Sequence[int],
            coords_origin: Sequence[float],
            spacing: float,
            data: Sequence[Any],
            data_labels: Optional[Sequence[Any]] = None,
            data_type: np.dtype = np.single,
    ):
        """
        Parameters
        ----------
        shape : Sequence[int]
            The number of grid points in each direction.
        data : Sequence[Any]
            Sequence of different data for each grid point. Each data should be in a form
            broadcastable to the grid with given `shape`, i.e. either a single value,
            or a sequence of values as many as the number of grid points. For sequential data,
            the values should be ordered such that the last index is moving fastest.
        data_labels : Sequence[str]
            Label for each type of data present. The grid can then be indexed using labels as well.
        coords_origin : Sequence[float]
            Real coordinates (e.g. in 3D this would be x, y, z) of the origin point,
            i.e. the point where all grid coordinates are zero.
            This is used to calculate the real coordinates of all other points.
        spacing : float
            The distance between two adjacent grid points. Since the grid is evenly spaced in
            all dimensions, this argument must be a single number.
            This is used to calculate the real coordinates of all other points.
        data_type : numpy.dtype, Optional, default: numpy.single
            The datatype of the grid point values
        """
        if len(coords_origin) != len(shape):
            raise ValueError("`coords_origin` and `shape` must have the same dimension.")
        # if data_labels is None:
        #     self._labels = np.arange(len(data))
        # elif len(data) == len(data_labels):
        #     self._labels = np.array(data_labels)
        # else:
        #     raise ValueError("`data` and `data_labels` must have the same length.")
        self._coords_origin = np.array(coords_origin)
        self._shape = np.array(shape)
        self._spacing = spacing

        self._tensor = np.empty(shape=(*self._shape, len(data)), dtype=data_type)
        for dat_idx, dat in enumerate(data):
            self._tensor[..., dat_idx] = dat
        return

    @property
    def shape(self):
        return tuple(self._shape)

    @property
    def spacing(self):
        return self._spacing

    # @property
    # def data_labels(self):
    #     return self._labels

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

    def __getitem__(self, item):
        return self._tensor.__getitem__(item)
        # def index_of_label(label) -> np.array:
        #     index = np.argwhere(self._labels == np.expand_dims(np.array(label), axis=-1))[:, -1]
        #     if index.size == 0:
        #         raise IndexError("data label not found.")
        #     elif index.size == 1:
        #         return index[0]
        #     return index
        # if isinstance(item, int) or (isinstance(item, tuple) and isinstance(item[-1], int)):
        #
        # if isinstance(item, str):
        #     return self._tensor.__getitem__(..., index_of_label(item))
        # elif isinstance(item, tuple) and isinstance(item[-1], (str, Sequence, np.ndarray)):
        #     return self._tensor.__getitem__(*item[:-1], index_of_label(item))
        # else:
