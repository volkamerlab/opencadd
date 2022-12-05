
# Standard library
from typing import Any, Literal, Sequence, Union, Optional, Tuple
import itertools
# 3rd-party
import numpy as np
import numpy.typing as npt
from opencadd.typing import ArrayLike
from opencadd.misc.spatial import Grid


class Field:
    """
    A collection of n scalar fields, or one n-dimensional vector field, sampled at regularly
    spaced points on a 3-dimensional grid in Euclidean space, over time (or e.g. in different
    environments).
    """

    @classmethod
    def from_array_like(
            cls,
            field_tensor: npt.ArrayLike,
            order_axes: Tuple[Literal[0, 1, 2, 3, 4]] = (0, 1, 2, 3, 4),
            field_datatype: npt.DTypeLike = np.single,
            grid_origin: Tuple[float, float, float] = (0, 0, 0),
            grid_point_spacing: float = 1,
    ):
        """
        Instantiate from a 5-dimensional array-like object containing the field values.

        Parameters
        ----------
        field_tensor : array_like
            A 5-dimensional array representing the field values for each grid point at different
            times. The dimensions represent time, index along x, y, and z axes, and field values.
            In each dimension, the elements should be ordered from the smallest index to largest.
        order_axes : Tuple[Literal[0, 1, 2, 3, 4]], Optional, default: (0, 1, 2, 3, 4)
            The axis index of time, x, y, and z axes, and field values, respectively.
        field_datatype : numpy.dtype, Optional, default: numpy.single
            The datatype of the field values.
        grid_origin : Tuple[float, float, float]
            Position (i.e. x, y, and z coordinates) of the origin grid point, i.e. the grid at
            index (0, 0, 0).
        grid_point_spacing : float
            The grid delta, i.e. distance between two adjacent grid points along x, y,
            and z axes. Since the grid is evenly spaced along all directions, this should be a
            scalar value. The unit should be the same as `position_origin`.

        Returns
        -------
        Field
        """
        if not isinstance(field_tensor, ArrayLike):
            raise ValueError("`field_values` must be array-like.")
        if not isinstance(order_axes, tuple) or len(order_axes) != 5:
            raise ValueError("`order_axes` must be array-like with length 5.")
        for axis in order_axes:
            if axis not in (0, 1, 2, 3, 4):
                raise ValueError("`order_axes` can only have integer values between 0 and 4.")
        # Create the underlying data structure and fill with data
        tensor: np.ndarray = np.array(field_tensor, dtype=field_datatype)
        if order_axes != (0, 1, 2, 3, 4):
            tensor = np.moveaxis(
                tensor,
                source=order_axes,
                destination=(0, 1, 2, 3, 4)
            )
        # len_field_values, len_field_names = len(field_values), len(field_names)
        # if len_field_values == 0:
        #     raise ValueError("`field_values` cannot be an empty array.")
        # if field_names is not None:
        #     if not isinstance(field_names, npt.ArrayLike):
        #         raise ValueError("`field_names` should either be None or array-like.")
        #     if len_field_names == 0:
        #         raise ValueError("`field_names` cannot be an empty sequence.")
        #     if len_field_names == 1:
        #         self._is_vector_field: bool = len_field_values > 1
        #     elif len_field_names != len_field_values:
        #         raise ValueError("`field_values` and `field_names` should have the same length.")
        #
        # # Set field type and name
        # self._field_names = field_names
        # if field_names is None:
        #     self._is_vector_field = True
        # elif len_field_names == 1:
        #     self._is_vector_field = len(field_values) > 1
        #
        # elif len(field_values) == len(field_names):
        #     self._field_names = np.array(field_names)
        # else:
        #     raise ValueError(
        #         "Number of provided field names does not match that of provided field values."
        #     )
        return cls(
            field_tensor=tensor,
            grid_origin=grid_origin,
            grid_point_spacing=grid_point_spacing
        )

    def __init__(
            self,
            field_tensor: np.ndarray,
            grid_origin: Tuple[float, float, float] = (0, 0, 0),
            grid_point_spacing: float = 1,
    ):
        """
        Parameters
        ----------
        field_tensor : array_like
            A 5-dimensional array representing the field values for each grid point at different
            times. The dimensions represent time, index along x, y, and z axes, and field values.
            In each dimension, the elements should be ordered from the smallest index to largest.

        field_names : Sequence[str]
            Label for each type of data present. The grid can then be indexed using labels as well.

        grid_origin : Tuple[float, float, float]
            Position (i.e. x, y, and z coordinates) of the origin grid point, i.e. the grid at
            index (0, 0, 0).
        grid_point_spacing : float
            The grid delta, i.e. distance between two adjacent grid points along x, y,
            and z axes. Since the grid is evenly spaced along all directions, this should be a
            scalar value. The unit should be the same as `position_origin`.
        """
        # Check for errors in input arguments.
        if not isinstance(field_tensor, np.ndarray):
            raise ValueError("`field_values` must be a numpy array.")
        if field_tensor.ndim != 5:
            raise ValueError("`field_values` must have 5 dimensions.")
        if not isinstance(grid_origin, ArrayLike) or len(grid_origin) != 3:
            raise ValueError("`position_origin` must be sequence of three integers.")
        if not isinstance(grid_origin, np.ndarray):
            for pos in grid_origin:
                if not isinstance(pos, (int, float)):
                    raise ValueError("`position_origin` must be sequence of three integers.")
        if not isinstance(grid_point_spacing, (float, int)) or grid_point_spacing <= 0:
            raise ValueError("`spacing` must be a positive number.")
        self._grid = Grid(
            shape=field_tensor.shape[1:4],
            origin=grid_origin,
            spacing=grid_point_spacing
        )
        self._tensor: np.ndarray = field_tensor
        return

    @property
    def grid(self) -> Grid:
        return self._grid

    @property
    def tensor(self):
        """
        The 5-dimensional array representing the field values for each grid point at different
        times.

        Returns
        -------
        numpy.ndarray
            A 5-dimensional array of shape (n_t, n_x, n_y, n_z, n_f), with
            n_t: temporal length of the field.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
            n_f: number of field values for each point.
        """
        return self._tensor

    def spatial_direction_vectors(self, dimensions=None):
        return np.pad(
            self._grid.direction_vectors(dimensions=dimensions),
            pad_width=((0, 0), (1, 0)),
            mode="constant",
            constant_values=0
        )

    @property
    def temporal_length(self) -> int:
        """
        Number of times the field values were sampled, e.g. at different times or different
        environments.
        """
        return self._tensor.shape[0]

    @property
    def fields_count(self) -> int:
        return self._tensor.shape[-1]


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


