from abc import ABC, abstractmethod
from typing import Sequence, Optional, Literal

import numpy as np
import jax.numpy as jnp

from opencadd._typing import ArrayLike
from opencadd.chemical import abc


class MolecularComposition (abc.Composition):
    def __init__(
            self,
            atomic_numbers: Sequence[int],
            mass_numbers: Optional[Sequence[int]] = None,
            charges: Optional[Sequence[float]] = None,
    ):
        self._z: jnp.ndarray
        self._a: jnp.ndarray
        self._e: jnp.ndarray

        self._z = jnp.asarray(atomic_numbers, dtype=np.ubyte)
        self._a = (
            jnp.asarray(mass_numbers, dtype=np.ubyte) if mass_numbers is not None
            else jnp.zeros_like(self._z)
        )
        self._e = (
            jnp.asarray(charges, dtype=np.single) if mass_numbers is not None
            else jnp.zeros_like(self._z)
        )
        return

    @property
    def atomic_numbers(self):
        return self._z

    @property
    def atomic_masses(self):
        return self._a

    @property
    def charges(self):
        return self._e

    def van_der_waals_radii(self, ref: Literal['default']):
        raise NotImplementedError


class ProteinComposition(MolecularComposition):
    pass