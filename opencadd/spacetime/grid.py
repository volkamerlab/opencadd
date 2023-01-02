"""
General functions and routines for spatial calculations, such as distances, volumes, etc.
"""


# Standard library
from typing import Any, Literal, Sequence, Union, Optional, Tuple
import itertools
# 3rd-party
import numpy as np
import jax.numpy as jnp
import scipy as sp
import numpy.typing as npt
from opencadd.typing import ArrayLike
from opencadd.spacetime.data_structure import TensorDataSet


def from_shape_spacing_center(
        center: Optional[Sequence[float]] = None,
):
    if origin is None:
        if center is None:
            self._origin = np.zeros(shape=self._dimension)
            self._center = self._size / 2
        else:
            self._center = np.array(center)
            self._origin = self._center - self._size / 2
    else:
        if center is None:
            self._origin = np.array(origin)
            self._center = self._origin + self._size / 2
        else:
            self._origin = np.array(origin)
            self._center = np.array(center)
            if not np.allclose(self._center, self._origin + self._size / 2):
                raise ValueError("Both `origin` and `center` are provided, but do not match.")
    pass


def from_shape_spacing_origin(
            shape: Sequence[int],
            spacing: float = 1,
            origin: Optional[Sequence[float]] = None,
):
    """
    Instantiate a grid by providing its shape, spacing, and optionally either origin or
    center coordinates. If neither is provided, the origin will be set to all zeros.

    Parameters
    ----------
    shape : Sequence[int]
        Shape of the grid, i.e. number of grid points per dimension.
    spacing : float, Optional, default: 1
        Distance between two adjacent points along a dimension.
    origin : Sequence[float], Optional, default: None
        Coordinates of the origin point of the grid, i.e. the point where all indices are zero.
    center : Sequence[float], Optional, default: None
        Coordinates of the geometric center of the grid.
    """



class Grid:
    """
    An n-dimensional grid of equidistant points in space.
    """

    __slots__ = (
        "_shape",
        "_spacing",
        "_meshgrid",
        "_tensor",
    )

    def __init__(self, ranges: npt.ArrayLike):
        """
        Parameters
        ----------
        ranges : ArrayLike, shape=(n, 3), dtype=[float, float, int]
            Start value, end value (inclusive), and total number of points each dimension.
        """
        # Declare attributes
        self._shape: jnp.ndarray
        self._spacing: jnp.ndarray
        self._meshgrid: jnp.ndarray
        self._tensor: TensorDataSet

        ranges_array = jnp.asarray(ranges)
        if ranges_array.ndim != 2 or ranges_array.shape[1] != 3:
            raise ValueError("Shape of `ranges` must be (n, 3).")
        if jnp.any(ranges_array[:,-1] <= 0):
            raise ValueError("Number of points in each dimension must be a positive integer.")

        self._shape = ranges_array[:, -1].astype(int)
        self._spacing = (
                jnp.absolute(ranges_array[:, 0] - ranges_array[:, 1])
                / (ranges_array[:, 2].astype(int) - 1)
        )

        # Calculating meshgrids with jax is sometimes inaccurate
        # (see https://github.com/google/jax/issues/13782).
        # Thus, create the meshgrid with numpy and transform to jax
        self._meshgrid = jnp.asarray(
            np.mgrid[
                tuple(map(lambda axis: slice(axis[0], axis[1], complex(0, int(axis[2]))), ranges))
            ]
        )

        self._tensor = TensorDataSet(
            data=jnp.expand_dims(jnp.stack(self._meshgrid, axis=-1), axis=0),
            title_instance="time",
            title_sample="grid",
            title_observation="position"
        )



        self._origin: np.ndarray
        self._center: np.ndarray
        self._dimension: int
        self._size: np.ndarray
        self._indices: np.ndarray
        self._shift_vectors: np.ndarray
        self._coordinates: np.ndarray
        self._direction_vectors: np.ndarray
        self._direction_vectors_dimension: np.ndarray

        # Assign attributes
        self._shape = np.array(shape).astype(int)
        self._spacing = spacing
        self._dimension = self._shape.size
        self._size = (self._shape - 1) * self._spacing
        self._indices = np.array(list(np.ndindex(*self._shape))).reshape(*self._shape, -1)

        if origin is None:
            if center is None:
                self._origin = np.zeros(shape=self._dimension)
                self._center = self._size / 2
            else:
                self._center = np.array(center)
                self._origin = self._center - self._size / 2
        else:
            if center is None:
                self._origin = np.array(origin)
                self._center = self._origin + self._size / 2
            else:
                self._origin = np.array(origin)
                self._center = np.array(center)
                if not np.allclose(self._center, self._origin + self._size / 2):
                    raise ValueError("Both `origin` and `center` are provided, but do not match.")

        self._shift_vectors = self._indices * self._spacing
        self._coordinates = self._shift_vectors + self._origin
        self._direction_vectors = np.array(
            list(itertools.product([-1, 0, 1], repeat=self._dimension))
        )
        self._direction_vectors_dimension = np.count_nonzero(self._direction_vectors, axis=-1)
        self._kdtree = sp.spatial.KDTree(data=self._coordinates.reshape(-1, self._dimension))
        return

    @property
    def dimension(self) -> int:
        """
        Dimension of the grid, i.e. number of axes.
        """
        return self._dimension

    @property
    def shape(self) -> Tuple[int]:
        """
        Shape of the grid, i.e. number of grid points in each dimension.
        """
        return tuple(self._shape)

    @property
    def size(self) -> Tuple[float]:
        """
        Size of the grid, i.e. length in each dimension,
        calculated from number of points and spacing.
        """
        return tuple(self._size)

    @property
    def point_count(self) -> int:
        """
        Total number of points in the grid.
        """
        return self._shape.prod()

    @property
    def origin(self) -> Tuple[float]:
        """
        Coordinates of the origin point of the grid, i.e. the point where all indices are zero.
        """
        return tuple(self._origin)

    @property
    def center(self) -> Tuple[float]:
        """
        Coordinates of the center of the grid.
        """
        return tuple(self._origin + self._size / 2)

    @property
    def spacing(self) -> float:
        """
        Distance between two adjacent points along a dimension.
        """
        return self._spacing

    @property
    def indices(self) -> np.ndarray:
        """
        Indices of all grid points.
        """
        return self._indices

    @property
    def coordinates(self) -> np.ndarray:
        """
        Coordinates of all grid points.

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
        return self._coordinates

    @property
    def coordinates_2d(self):
        return

    @property
    def kdtree(self) -> sp.spatial.KDTree:
        return self._kdtree

    def direction_vectors(self, dimensions: Optional[Sequence[int]] = None) -> np.ndarray:
        if dimensions is None:
            dimensions = np.arange(1, self._dimension + 1)
        return self._direction_vectors[np.isin(self._direction_vectors_dimension, dimensions)]

    def distance(
            self,
            coordinates: np.ndarray
    ):
        """
        Calculate distances between each grid point and a set of points.

        Parameters
        ----------
        coordinates : numpy.ndarray
            A 2-dimensional array of shape (n, d), containing the coordinates of n points in a
            d-dimensional space, where d is equal to the grid dimension.

        Returns
        -------
        dists : numpy.ndarray
            Distances between each grid point and each point in `coordinates`.
        """

        dist_vects_origin = coordinates - self._origin
        dist_vects_origin_repeated = np.tile(
            dist_vects_origin[np.newaxis],
            reps=(self.point_count, 1, 1) if coordinates.ndim == 2 else (1, self.point_count, 1, 1)
        ).reshape(
            (*self.shape, *coordinates.shape) if coordinates.ndim == 2
            else (coordinates.shape[0], *self.shape, *coordinates.shape[1:])
        )
        dist_vects = dist_vects_origin_repeated - self._shift_vectors[..., np.newaxis, :]
        return np.linalg.norm(dist_vects, axis=-1)




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


def xeno_neighbor_distance(
        bool_array: np.ndarray,
        dir_vectors: np.ndarray,
        dir_multipliers: Optional[Union[Sequence[int], int]] = None,
):
    """
    Given an n-dimensional boolean array, for each boolean element calculate its distances
    to the first opposite elements along a number of given directions.

    Parameters
    ----------
    bool_array : numpy.ndarray
        An n-dimensional boolean array.
    dir_vectors : numpy.ndarray
        A 2-dimensional array of shape (k, n), containing k direction vectors in an
        n-dimensional space.
    dir_multipliers : Sequence[int] | int, Optional, default: None
        Maximum multipliers for direction vectors, i.e. maximum number of times to travel along
        each direction to find the opposite neighbor of an element, before terminating the
        search. This can be a single integer used for all direction vectors, or a sequence of k
        integers, one for each direction vector. If not provided, search will continue until one
        edge of the array is reached.

    Returns
    -------
    numpy.ndarray
        An (n+1)-dimensional array of integers, where the first n dimensions match the shape of
        the input `bool_array`, and the last dimension has k elements, each describing the
        distance to the nearest opposite element along the corresponding direction vector in
        `dir_vectors`. The values are all integers, and correspond to the number of times to
        travel along the corresponding direction vector to reach the nearest xeno element.
        The value will be 0 for directions where no opposite element was found.
    """
    def slicer(vec):
        """
        For a given displacement vector (i.e. direction vector times a multiplier), calculate
        two tuples of slices, which index the starting elements and end elements of that
        displacement on the array.
        """
        start_slices, end_slices = [], []
        for val in vec:
            if val > 0:
                start_slices.append(slice(None, -val))
                end_slices.append(slice(val, None))
            elif val < 0:
                start_slices.append(slice(-val, None))
                end_slices.append(slice(None, val))
            else:
                for lis in [start_slices, end_slices]:
                    lis.append(slice(None, None))
        return tuple(start_slices), tuple(end_slices)
    # Initiate the array of distance with zeros.
    dists = np.zeros(shape=(*bool_array.shape, dir_vectors.shape[0]), dtype=np.uintc)
    # Calculate the maximum multiplier along each direction:
    # First, calculate the maximum possible multipliers
    with np.errstate(divide='ignore', invalid='ignore'):
        max_mult_axis = (np.array(bool_array.shape) - 1) / np.abs(dir_vectors)
        max_mult_axis[np.isnan(max_mult_axis)] = np.inf
        max_mult_dir = np.min(max_mult_axis, axis=-1)
    # Then, compare with user-input multipliers and take the smaller one in each direction.
    if dir_multipliers is None:
        dir_multipliers = np.ones(dir_vectors.shape[0]) * np.max(bool_array.shape)
    elif isinstance(dir_multipliers, int):
        dir_multipliers = np.ones(dir_vectors.shape[0]) * dir_multipliers
    max_mult = np.min((max_mult_dir, dir_multipliers), axis=0).astype(int) + 1
    # Loop through directions, and for each direction through multipliers, and calculate
    # distances between starting elements and end elements.
    for idx_dir, direction in enumerate(dir_vectors):
        curr_mask = np.ones_like(bool_array)
        for mult in range(1, max_mult[idx_dir]):
            start_slice, end_slice = slicer(mult*direction)
            reached_xeno = np.logical_xor(bool_array[start_slice], bool_array[end_slice])
            dists[(*start_slice, idx_dir)][curr_mask[start_slice]] = reached_xeno[curr_mask[
                start_slice]] * mult
            curr_mask[start_slice][reached_xeno] = 0
    return dists


def from_bounds_and_spacing(
        lower_bounds: npt.ArrayLike,
        upper_bounds: npt.ArrayLike,
        spacings: npt.ArrayLike,
) -> Grid:
    return