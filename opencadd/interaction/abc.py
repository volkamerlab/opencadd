"""
Abstract base clases used in the package.
"""

from abc import abstractmethod, ABC
from typing import Tuple
import numpy as np
from opencadd.misc import field


class IntraMolecularInteractionField(ABC, field.Field):
    """
    Intramolecular interaction field.
    """
    def __init__(
            self,
            field_tensor: np.ndarray,
            grid_origin: Tuple[float, float, float] = (0, 0, 0),
            grid_point_spacing: float = 1,
    ):
        super().__init__(
            field_tensor=field_tensor,
            grid_origin=grid_origin,
            grid_point_spacing=grid_point_spacing
        )
        return

    @property
    @abstractmethod
    def electrostatic(self):
        ...

    @property
    @abstractmethod
    def van_der_waals(self):
        ...

    @property
    @abstractmethod
    def h_bond_donor(self):
        ...

    @property
    @abstractmethod
    def h_bond_acceptor(self):
        ...

    @property
    @abstractmethod
    def hydrophobic(self):
        ...

    @property
    @abstractmethod
    def aromatic(self):
        ...
