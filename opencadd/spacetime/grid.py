"""
An n-dimensional grid of points in euclidean space.
"""


# Standard library
from typing import Any, Literal, Sequence, Union, Optional, Tuple
import itertools
# 3rd-party
import numpy as np
import jax.numpy as jnp
import scipy as sp
import numpy.typing as npt
from opencadd._typing import ArrayLike
import opencadd as oc


class Grid:
    """
    An n-dimensional grid of points in euclidean space.
    """

    def __init__(
            self,
            shape: np.ndarray,
            size: np.ndarray,
            lower_bounds: np.ndarray,
            center: np.ndarray,
            upper_bounds: np.ndarray,
            spacings: np.ndarray,
            mgrid: np.ndarray,
    ):
        """

        Parameters
        ----------
        shape : numpy.ndarray
            Shape of the grid, i.e. number of points in each dimension.
        size : numpy.ndarray
            Length of the grid in each dimension.
        lower_bounds : numpy.ndarray
            Coordinates of the point with minimum values in all dimensions.
        center : numpy.ndarray
            Coordinates of the geometric center of the grid.
        upper_bounds : numpy.ndarray
            Coordinates of the point with maximum values in all dimensions.
        spacings : numpy.ndarray
            Spacing between grid points in each dimension.
        mgrid : numpy.ndarray
            Fleshed out meshgrid of grid point coordinates.
        """
        self._shape: np.ndarray = shape
        self._size: np.ndarray = size
        self._lower_bounds: np.ndarray = lower_bounds
        self._center: np.ndarray = center
        self._upper_bounds: np.ndarray = upper_bounds
        self._spacings: np.ndarray = spacings
        self._mgrid: jnp.ndarray = jnp.asarray(mgrid)

        self._dimension: int = self._shape.size
        self._coordinates: jnp.ndarray = jnp.stack(mgrid, axis=-1)
        self._count_points = np.prod(self._shape)
        self._indices: np.ndarray = np.array(list(np.ndindex(*self._shape))).reshape(*self._shape, -1)
        self._pointcloud = oc.spacetime.pointcloud.DynamicPointCloud(
            data=self._coordinates.reshape(1, self._count_points, self._dimension)
        )
        self._direction_vectors = np.array(
            list(itertools.product([-1, 0, 1], repeat=self._dimension))
        )
        self._direction_vectors_dimension = np.count_nonzero(self._direction_vectors, axis=-1)
        return

    @property
    def dimension(self) -> int:
        """
        Dimension of the grid, i.e. number of axes.
        """
        return self._dimension

    @property
    def shape(self) -> np.ndarray:
        """
        Shape of the grid, i.e. number of grid points in each dimension.
        """
        return np.array(self._shape)

    @property
    def size(self) -> np.ndarray:
        """
        Size of the grid, i.e. length in each dimension.
        """
        return np.array(self._size)

    @property
    def count_points(self) -> int:
        """
        Total number of points in the grid.
        """
        return self._count_points

    @property
    def lower_bounds(self) -> np.ndarray:
        """
        Coordinates of the origin point of the grid, i.e. the point where all indices are zero.
        """
        return np.array(self._lower_bounds)

    @property
    def center(self) -> np.ndarray:
        """
        Coordinates of the center of the grid.
        """
        return np.array(self._center)

    @property
    def upper_bounds(self) -> np.ndarray:
        return np.array(self._upper_bounds)

    @property
    def spacings(self) -> np.ndarray:
        """
        Distance between two adjacent points along a dimension.
        """
        return np.array(self._spacings)

    @property
    def indices(self) -> np.ndarray:
        """
        Indices of all grid points.
        """
        return self._indices

    @property
    def coordinates(self) -> jnp.ndarray:
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
    def points(self):
        return self._pointcloud

    def direction_vectors(self, dimensions: Optional[Sequence[int]] = None) -> np.ndarray:
        if dimensions is None:
            dimensions = np.arange(1, self._dimension + 1)
        return self._direction_vectors[np.isin(self._direction_vectors_dimension, dimensions)]


def from_bounds_shape(
        lower_bounds: Sequence[float],
        upper_bounds: Sequence[float],
        shape: Sequence[int],
) -> Grid:
    """
    Create a `Grid` from its lower- and upper bounds, and shape.

    Parameters
    ----------
    lower_bounds : sequence of float
        Coordinates of the point with minimum values in all dimensions.
    upper_bounds : sequence of float
        Coordinates of the point with maximum values in all dimensions.
        This must have the same length as `lower_bounds`.
    shape : sequence of int
        Shape of the grid, i.e. number of points in each dimension.
        This must have the same length as `lower_bounds` and `upper_bounds`.

    Returns
    -------
    Grid
    """
    lower_bounds = np.asarray(lower_bounds)
    upper_bounds = np.asarray(upper_bounds)
    shape = np.asarray(shape)
    for bound, arg_name in zip((lower_bounds, upper_bounds, shape), ("lower_bounds", "upper_bounds", "shape")):
        if bound.ndim != 1:
            raise ValueError(
                f"Parameter `{arg_name}` expects a 1D array, "
                f"but input argument had {bound.ndim} dimensions. Input was: {bound}"
            )
    for bound, arg_name in zip((lower_bounds, upper_bounds), ("lower_bounds", "upper_bounds")):
        if not (np.issubdtype(bound.dtype, np.floating) or np.issubdtype(bound.dtype, np.integer)):
            raise ValueError(
                f"Parameter `{arg_name}` expects an array of real numbers, "
                f"but input argument had elements of type {bound.dtype}. Input was: {bound}"
            )
    if not np.issubdtype(shape.dtype, np.integer):
        raise ValueError(
            f"Parameter `shape` expects an array of integers, "
            f"but input argument had elements of type {shape.dtype}. Input was: {shape}"
        )
    for arg, arg_name in zip((upper_bounds, shape), ("upper_bounds", "shape")):
        if lower_bounds.size != arg.size:
            raise ValueError(
                f"Parameters `lower_bound` and `{arg_name}` expect 1D arrays of same size, "
                f"but input argument `lower_bounds` had size {lower_bounds.size}, "
                f"while `{arg_name}` had size {arg.size}. "
                f"Inputs were: `lower_bounds` = {lower_bounds}\n`{arg_name}` = {arg}."
            )
    size = upper_bounds - lower_bounds
    size_is_invalid = size <= 0
    if np.any(size_is_invalid):
        raise ValueError(
            "All values in `lower_bounds` must be strictly smaller "
            f"than corresponding values in `upper_bounds`, but at indices {np.where(size_is_invalid)[0]} "
            f"`lower_bounds` had values {lower_bounds[size_is_invalid]} and `upper_bounds` had "
            f"values {upper_bounds[size_is_invalid]}."
        )
    shape_is_invalid = shape <= 0
    if np.any(shape_is_invalid):
        raise ValueError(
            "All values in `shape` must be positive integers, "
            f"but at indices {np.where(shape_is_invalid)[0]} "
            f"`shape` had values {shape[shape_is_invalid]}."
        )
    slices = tuple(
        slice(start, end, complex(num_points))
        for start, end, num_points in zip(lower_bounds, upper_bounds, shape)
    )
    return Grid(
        shape=shape,
        size=size,
        lower_bounds=lower_bounds,
        center=(lower_bounds + upper_bounds) / 2,
        upper_bounds=upper_bounds,
        spacings=size / (shape - 1),
        mgrid=np.mgrid[slices]
    )


def from_bounds_spacing(
        lower_bounds: Sequence[float],
        upper_bounds: Sequence[float],
        spacings: Sequence[float],
) -> Grid:
    """
    Create a `Grid` from its lower- and upper bounds, and spacings.

    Parameters
    ----------
    lower_bounds : sequence of float
        Coordinates of the point with minimum values in all dimensions.
    upper_bounds : sequence of float
        Coordinates of the point with maximum values in all dimensions.
        This must have the same length as `lower_bounds`.
    spacings : sequence of int
        Spacing between grid points in each dimension.
        This must have the same length as `lower_bounds` and `upper_bounds`.

    Returns
    -------
    Grid
    """
    lower_bounds = np.asarray(lower_bounds)
    upper_bounds = np.asarray(upper_bounds)
    spacings = np.asarray(spacings)
    size = upper_bounds - lower_bounds
    shape = (size / spacings) + 1
    return from_bounds_shape(lower_bounds=lower_bounds, upper_bounds=upper_bounds, shape=shape)


def from_center_spacing_size(
        center: Sequence[float],
        spacings: Sequence[float],
        size: Sequence[float],
) -> Grid:
    """
    Create a `Grid` from its center, size, and spacings.

    Parameters
    ----------
    center : sequence of float
        Coordinates of the geometric center of the grid.
    spacings : sequence of int
        Spacing between grid points in each dimension.
        This must have the same length as `center` and `size`.
    size : sequence of float
        Length of the grid in each dimension.
        This must have the same length as `center`.

    Returns
    -------
    Grid
    """
    center = np.asarray(center)
    size = np.asarray(size)
    spacings = np.asarray(spacings)
    return from_bounds_spacing(
        lower_bounds=center-size/2,
        upper_bounds=center+size/2,
        spacings=spacings
    )


def from_center_spacing_shape(
        center: Sequence[float],
        spacings: Sequence[float],
        shape: Sequence[int],
):
    center = np.asarray(center)
    shape = np.asarray(shape)
    spacings = np.asarray(spacings)
    return from_center_spacing_size(center=center, spacings=spacings, size=(shape - 1)*spacings)


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

