from abc import ABC, abstractmethod
from typing import Sequence
import ipywidgets

from opencadd.spacetime.abc import Volume


class BindingSite(Volume, ABC):

    @abstractmethod
    @property
    def viewer(self) -> ipywidgets.Widget:
        """
        An interactive Jupyter widget for viewing the binding site.

        Returns
        -------
        ipywidgets.Widget
        """
        ...

