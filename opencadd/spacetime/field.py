
# Standard library
from typing import Literal, Sequence, Optional, Tuple, Union, Any
import operator
# 3rd-party
import numpy as np
import numpy.typing as npt
import jax.numpy as jnp
from opencadd._typing import ArrayLike
from opencadd.spacetime.grid import Grid
from opencadd import _exceptions


class ToxelField:
    """
    A collection of n scalar fields, or one n-dimensional vector field, sampled at regularly
    spaced points on a 3-dimensional grid in Euclidean space, over time (or e.g. in different
    environments).
    """

    def __init__(
            self,
            tensor: ArrayLike,
            grid: Grid,
            names: Optional[Sequence[Any]] = None,
    ):
        """
        Parameters
        ----------
        tensor : array_like
            A 5-dimensional array representing the field values for each grid point at different
            times. The dimensions represent time, index along x, y, and z axes, and field values.
            In each dimension, the elements should be ordered from the smallest index to largest.

        names : Sequence[str]
            Label for each type of data present. The grid can then be indexed using labels as well.
        """
        # Check for errors in input arguments.
        _exceptions.raise_type(self.__class__.__name__, ("grid", grid, Grid))
        tensor = jnp.asarray(tensor)
        _exceptions.raise_array(
            parent_name=self.__class__.__name__,
            param_name="tensor",
            array=tensor,
            ndim_gt=2
        )
        if grid.dimension != tensor.ndim - 2:
            raise ValueError(
                "Dimension Mismatch:\n"
                "An n-dimensional `tensor` requires an (n - 2)-dimensional `grid`, "
                f"but the input tensor was {tensor.ndim}-dimensional, while the grid "
                f"was {grid.dimension}-dimensional."
            )
        if np.any(grid.shape != tensor.shape[1:-1]):
            raise ValueError(
                "Shape Mismatch:\n"
                "The spatial shape of the input tensor (i.e. `tensor.shape[1:-1]`) "
                f"must be equal to the shape of the grid, but the input tensor had "
                f"a spatial shape of {tensor.shape[1:-1]}, while the shape of the grid "
                f"was {grid.shape}."
            )
        if names is None:
            names = np.arange(tensor.shape[-1])
        else:
            names = np.asarray(names)
            _exceptions.raise_array(
                parent_name=self.__class__.__name__,
                param_name="names",
                array=names,
                ndim_eq=1
            )
            if names.size != tensor.shape[-1]:
                raise ValueError(
                    "Shape Mismatch:\n"
                    "The size of `names` must be equal to the shape of `tensor` along its last dimension, "
                    f"but the input array had {names.size} elements, while the last dimension of `tensor` "
                    f"had {tensor.shape[-1]} elements."
                )

        self._tensor: jnp.ndarray = jnp.asarray(tensor)
        self._grid: Grid = grid
        self._field_names = names
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

    def __call__(self, name=None, instance=slice(None)):
        idx_fields = slice(None) if name is None else self.index_field(name)
        return self._tensor[instance, ..., idx_fields]

    @property
    def temporal_length(self) -> int:
        """
        Number of times the field values were sampled, e.g. at different times or different
        environments.
        """
        return self._tensor.shape[0]

    @property
    def field_names(self):
        return self._field_names

    @property
    def fields_count(self) -> int:
        return self._tensor.shape[-1]

    def index_field(self, name: Union[Any, Sequence[Any]]) -> np.ndarray:
        names = np.asarray(name).reshape(-1, 1)
        indices = np.argwhere(self._field_names == names)
        if indices.shape[0] != names.size:
            ind_bad_names = np.setdiff1d(np.arange(names.size), indices[:, 0])
            raise IndexError(
                f"Following field names are not valid: {names[ind_bad_names]}. "
                f"Valid field names are: {self.field_names}."
            )
        return np.squeeze(indices[:, 1])


    def calculate_vacancy(
            self,
            energy_cutoff: float = +0.6,
            mode: Optional[Literal["max", "min", "avg", "sum"]] = "min",
    ) -> np.ndarray:
        """
        Calculate whether each grid point is vacant, or occupied by a target atom.

        Parameters
        ----------
        energy_cutoff : float, Optional, default: +0.6
            Cutoff value for energy; grid points with energies lower than cutoff are considered
            vacant.
        mode: Literal["max", "min", "avg", "sum"], Optional, default: "min"
            If the energy of more than one ligand type is to be compared, this parameter defines
            how those different energy values must be processed, before comparing with the cutoff.
        ligand_types : Sequence[opencadd.consts.autodock.AtomType], Optional, default: None
            A subset of ligand types that were used to initialize the object, whose energy values
            are to be taken as reference for calculating the vacancy of each grid point. If not
            set to None, then all ligand interaction energies are considered.

        Returns
        -------
        vacancy : numpy.ndarray[dtype=numpy.bool_, shape=T2FPharm.grid.shape[:-1]]
            A 4-dimensional boolean array matching the first four dimensions of `T2FPharm.grid`,
            indicating whether each grid point is vacant (True), or occupied (False).
            Vacant grid points can easily be indexed by `T2FPharm.grid[vacancy]`.
        """
        # The reducing operations corresponding to each `mode`:
        red_fun = {"max": np.max, "min": np.min, "avg": np.mean, "sum": np.sum}
        # Get index of input ligand types
        # if ligand_types is None:
        #     ind = slice(None)
        # else:
        #     ind = np.argwhere(np.expand_dims(ligand_types, axis=-1) == self._probe_types)[:, 1]
        #     # Verify that all input ligand types are valid
        #     if len(ind) != len(ligand_types):
        #         raise ValueError(f"Some of input energies were not calculated.")
        # Reduce the given references using the given operation.
        energy_vals = red_fun[mode](self._interaction_field.van_der_waals, axis=-1)
        # Apply cutoff and return
        self._vacancy = energy_vals < energy_cutoff
        return self._vacancy

    def __getitem__(self, item):
        return self._tensor.__getitem__(item)
        # if isinstance(item, int) or (isinstance(item, tuple) and isinstance(item[-1], int)):
        #
        # if isinstance(item, str):
        #     return self._tensor.__getitem__(..., index_of_label(item))
        # elif isinstance(item, tuple) and isinstance(item[-1], (str, Sequence, np.ndarray)):
        #     return self._tensor.__getitem__(*item[:-1], index_of_label(item))
        # else:

    def visualize(self):
        pass


def from_tensor_grid(
        tensor,
        grid,
        names,
):
    return ToxelField(tensor=tensor, grid=grid, names=names)


def from_array_like(
        field_tensor: npt.ArrayLike,
        order_axes: Tuple[Literal[0, 1, 2, 3, 4]] = (0, 1, 2, 3, 4),
        field_names: Optional[Sequence[str]] = None,
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

    grid_origin : Tuple[float, float, float]
            Position (i.e. x, y, and z coordinates) of the origin grid point, i.e. the grid at
            index (0, 0, 0).
    grid_point_spacing : float
        The grid delta, i.e. distance between two adjacent grid points along x, y,
        and z axes. Since the grid is evenly spaced along all directions, this should be a
        scalar value. The unit should be the same as `position_origin`.

    Returns
    -------
    ToxelField
    """
    if not isinstance(field_tensor, ArrayLike):
        raise ValueError("`field_values` must be array-like.")
    if not isinstance(order_axes, tuple) or len(order_axes) != 5:
        raise ValueError("`order_axes` must be array-like with length 5.")
    for axis in order_axes:
        if axis not in (0, 1, 2, 3, 4):
            raise ValueError("`order_axes` can only have integer values between 0 and 4.")

    if not isinstance(grid_origin, ArrayLike) or len(grid_origin) != 3:
        raise ValueError("`position_origin` must be sequence of three integers.")
    if not isinstance(grid_origin, np.ndarray):
        for pos in grid_origin:
            if not isinstance(pos, (int, float)):
                raise ValueError("`position_origin` must be sequence of three integers.")
    if not isinstance(grid_point_spacing, (float, int)) or grid_point_spacing <= 0:
        raise ValueError("`spacing` must be a positive number.")


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
        field_names=field_names,
        grid_origin=grid_origin,
        grid_point_spacing=grid_point_spacing
    )

def empty(
        self,
        temporal_length: int,
        grid_shape: Tuple[int, int, int],
        grid_spacing: float,
        field_shape: Tuple[int],
        field_datatype: npt.DTypeLike = np.single,
        field_names: Optional[Sequence[str]] = None,
        grid_origin: Tuple[float, float, float] = (0, 0, 0),
):
    pass