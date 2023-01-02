from typing import Sequence
from abc import ABC, abstractmethod

import numpy as np
import jax.numpy as jnp
import numpy.typing as npt


class Volume(ABC):
    """
    An n-dimensional volume form (i.e. length, area, volume, hyper-volume), sampled at one or
    several instances.
    """

    @property
    @abstractmethod
    def size(self) -> npt.ArrayLike:
        """
        Size of the volume at each instance.
        """
        ...

    @property
    @abstractmethod
    def aabb_lower_bounds(self) -> npt.ArrayLike:
        """
        Lower bounds of the axis-aligned minimum bounding box.

        Returns
        -------
        ArrayLike,
        """
        ...

    @property
    @abstractmethod
    def aabb_upper_bounds(self) -> npt.ArrayLike:
        """
        Upper bounds of the axis-aligned minimum bounding box.
        """

    @property
    @abstractmethod
    def aabb_lengths(self) -> npt.ArrayLike:
        """
        Length
        """
