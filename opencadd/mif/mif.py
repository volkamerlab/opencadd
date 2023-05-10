"""
Abstract base clases used in the package.
"""

from typing import Tuple, Sequence
import numpy as np
from opencadd import spacetime
from opencadd import chem


class MolecularInteractionField:
    """
    Intramolecular interaction field.
    """
    def __init__(
            self,
            ensemble: chem.ensemble.ChemicalEnsemble,
            field: spacetime.field.ToxelField,
    ):
        self._ensemble = ensemble
        self._field = field
        return

    @property
    def ensemble(self):
        return self._ensemble

    @property
    def field(self):
        return self._field
