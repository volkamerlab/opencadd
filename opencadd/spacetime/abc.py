from typing import Sequence
from abc import ABC, abstractmethod


class Volume(ABC):
    """
    An n-dimensional volume form (i.e. length, area, volume, hyper-volume), sampled at one or
    several instances.
    """

    @property
    @abstractmethod
    def size(self) -> Sequence[float]:
        """
        Size of the volume at each instance.
        """
        ...

    @property
    @abstractmethod
    def minimum_bounding_box(self):
        ...


class Field(ABC):
    pass


class Grid(ABC):
    pass




