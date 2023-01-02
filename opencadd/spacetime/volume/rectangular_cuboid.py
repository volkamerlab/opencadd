import opencadd as oc
from opencadd.spacetime.volume import abc

import numpy as np
import numpy.typing as npt
import jax.numpy as jnp


class RectangularCuboid(abc.Volume):
    """
    An n-dimensional rectangular cuboid (i.e. line, rectangle, rectangular cuboid,
    hyper-rectangular cuboid), sampled at one or several instances.
    """

    def __init__(
            self,
            lower_bounds: npt.ArrayLike,
            upper_bounds: npt.ArrayLike,
    ):
        # Declare attributes
        self._lower_bounds: jnp.ndarray
        self._upper_bounds: jnp.ndarray
        # Check for errors in input
        lower_bounds_array = jnp.asarray(lower_bounds)
        upper_bounds_array = jnp.asarray(upper_bounds)
        if lower_bounds_array.ndim != 2 or upper_bounds_array != 2:
            raise ValueError("Input arrays must be 2-dimensional.")
        if lower_bounds_array.shape != upper_bounds_array.shape:
            raise ValueError("Input arrays must have the same shape.")
        # Assign attributes
        self._lower_bounds = lower_bounds_array
        self._upper_bounds = upper_bounds_array
        return

    @property
    def lower_bounds(self) -> jnp.ndarray:
        return self._lower_bounds

    @property
    def upper_bounds(self) -> jnp.ndarray:
        return self._upper_bounds



