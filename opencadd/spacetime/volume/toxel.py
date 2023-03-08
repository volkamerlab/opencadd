from typing import Optional, Union, Sequence

import numpy as np
import jax.numpy as jnp


from opencadd.spacetime.volume import abc
import opencadd as oc
import opencadd._typing
from opencadd import _exceptions


class ToxelVolume:

    def __init__(self, toxels, grid):
        _exceptions.raise_type(self.__class__.__name__, ("grid", grid, oc.spacetime.grid.Grid))
        toxels = jnp.asarray(toxels)
        _exceptions.raise_array(
            parent_name=self.__class__.__name__,
            param_name="toxels",
            array=toxels,
            ndim_gt=2,
            dtype=np.bool_,
        )
        if grid.dimension != toxels.ndim - 1:
            raise ValueError(
                "Dimension Mismatch:\n"
                "An n-dimensional `toxels` array requires an (n - 1)-dimensional `grid`, "
                f"but the input array was {toxels.ndim}-dimensional, while the grid "
                f"was {grid.dimension}-dimensional."
            )
        if np.any(grid.shape != toxels.shape[1:]):
            raise ValueError(
                "Shape Mismatch:\n"
                "The spatial shape of the toxels (i.e. `toxels.shape[1:]`) "
                f"must be equal to the shape of the grid, but the input array had "
                f"a spatial shape of {toxels.shape[1:]}, while the shape of the grid "
                f"was {grid.shape}."
            )
        self._toxels = jnp.asarray(toxels)
        self._grid = grid
        return

    @property
    def grid(self):
        return self._grid

    @property
    def toxels(self):
        return self._toxels

    def xeno_neighbor_distance(
            self,
            dir_vectors: Optional[np.ndarray] = None,
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

        if dir_vectors is None:
            dir_vectors = self.spatial_direction_vectors()
        # Initiate the array of distance with zeros.
        dists = np.zeros(shape=(*self._toxels.shape, dir_vectors.shape[0]), dtype=np.uintc)
        # Calculate the maximum multiplier along each direction:
        # First, calculate the maximum possible multipliers
        with np.errstate(divide='ignore', invalid='ignore'):
            max_mult_axis = (np.array(self._toxels.shape) - 1) / np.abs(dir_vectors)
            max_mult_axis[np.isnan(max_mult_axis)] = np.inf
            max_mult_dir = np.min(max_mult_axis, axis=-1)
        # Then, compare with user-input multipliers and take the smaller one in each direction.
        if dir_multipliers is None:
            dir_multipliers = np.ones(dir_vectors.shape[0]) * np.max(self._toxels.shape)
        elif isinstance(dir_multipliers, int):
            dir_multipliers = np.ones(dir_vectors.shape[0]) * dir_multipliers
        max_mult = np.min((max_mult_dir, dir_multipliers), axis=0).astype(int) + 1
        # Loop through directions, and for each direction through multipliers, and calculate
        # distances between starting elements and end elements.
        for idx_dir, direction in enumerate(dir_vectors):
            curr_mask = np.ones_like(self._toxels)
            for mult in range(1, max_mult[idx_dir]):
                start_slice, end_slice = slicer(mult * direction)
                reached_xeno = np.logical_xor(self._toxels[start_slice], self._toxels[end_slice])
                dists[(*start_slice, idx_dir)][curr_mask[start_slice]] = reached_xeno[curr_mask[
                    start_slice]] * mult
                curr_mask[start_slice][reached_xeno] = 0
        return dists

    def spatial_direction_vectors(self, dimensions=None):
        return np.pad(
            self._grid.direction_vectors(dimensions=dimensions),
            pad_width=((0, 0), (1, 0)),
            mode="constant",
            constant_values=0
        )