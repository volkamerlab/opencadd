from typing import Sequence, Union, List

import numpy as np
import scipy as sp
import jax.numpy as jnp


from opencadd._typing import ArrayLike
from opencadd.spacetime.spatial import Spatial


class Trajectory:
    def __init__(self, coordinates: ArrayLike):
        self._coordinates: jnp.ndarray
        self._spatials: List[Spatial]

        coordinates_array = jnp.asarray(coordinates)
        return
