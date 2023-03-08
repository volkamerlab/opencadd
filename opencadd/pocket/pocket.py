from typing import Sequence

import numpy as np
from opencadd.spacetime.volume import ToxelVolume


class BindingPocket:
    def __init__(
            self,
            ensemble,
            volume: ToxelVolume,
            atoms: Sequence[int],
    ):
        self._ensemble = ensemble
        self._volume = volume
        self._atom_ids = np.asarray(atoms)
        return

    @property
    def ensemble(self):
        return self._ensemble

    @property
    def volume(self):
        return self._volume

    def viewer(self):
        """
        An interactive Jupyter widget for viewing the binding site.

        Returns
        -------
        ipywidgets.Widget
        """
        return