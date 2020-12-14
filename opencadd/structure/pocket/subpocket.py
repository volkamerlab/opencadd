"""
opencadd.structure.pocket.subpocket

Defines subpockets.
"""

import logging

import numpy as np
import pandas as pd

_logger = logging.getLogger(__name__)


class Subpocket:
    """
    Class defining a subpocket.

    Attributes
    ----------
    name : str
        Subpocket name.
    color : str
        Region color name (matplotlib name).
    center : np.array
        Coordinates (x, y, z) of the subpocket center,
        i.e. the centroid of all anchor residues' CA atoms.
    _anchor_residues : list of Residue
        List of anchor residues.
    """

    def __init__(self, anchor_residues, name, color="blue", pocket_name=None):
        """
        Initialize a Subpocket object.

        Parameters
        ----------
        anchor_residues : list of Residue
            List of anchor residues.
        name : str
            Subpocket name.
        color : str
            Subpocket color (matplotlib name), blue by default.
        _pocket_name : str or None
            Pocket name.
        """

        self._anchor_residues = anchor_residues
        self.name = name
        self.color = color
        self._pocket_name = pocket_name
        self.center = self._centroid()

    @property
    def subpocket(self):
        """
        Subpocket attributes.

        Returns
        -------
        pd.Series
            Subpocket attributes.
        """
        return pd.Series(
            {
                "subpocket.name": self.name,
                "subpocket.color": self.color,
                "subpocket.center": self.center,
            }
        )

    @property
    def anchor_residues(self):
        """
        Anchor residues for all subpockets.
        - Subpocket name and color
        - Anchor residue IDs (user-defined input IDs or alternative
          IDs if input was not available)
        - Anchor residue indices
        - Anchor residue centers (coordinates)

        Returns
        -------
        pd.DataFrame
            Anchor residues for all subpockets.
        """

        anchor_residues = pd.DataFrame(
            [anchor_residue.anchor_residue for anchor_residue in self._anchor_residues]
        )
        anchor_residues.insert(0, "subpocket.name", self.name)
        return anchor_residues

    def _centroid(self):
        """
        Calculate the centroid of given input anchor residue centers.

        Returns
        -------
        np.array
            Subpocket center, i.e. the centroid of all anchor residue centers.
            None if anchor residues are missing.
        """

        anchor_residue_centers = [
            anchor_residue.center for anchor_residue in self._anchor_residues
        ]
        # Are there empty anchor residue centers?
        anchor_residue_centers_none = [
            center for center in anchor_residue_centers if center is None
        ]
        # If so, do not return a subpocket center.
        if len(anchor_residue_centers_none) != 0:
            return None
        # Else, calculate the centroid of all given anchor residue centers.
        else:
            subpocket_center = np.mean(anchor_residue_centers, axis=0)
            return subpocket_center
