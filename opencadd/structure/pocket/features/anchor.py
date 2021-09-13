"""
opencadd.structure.pocket.anchor

Defines anchor residues.
"""

import logging

import numpy as np
import pandas as pd

_logger = logging.getLogger(__name__)


class AnchorResidue:
    """
    Class defining a subpocket anchor.

    Attributes
    ----------
    center : numpy.array
        Coordinates (x, y, z) of the residue center.
    residue_id : str
        Residue ID.
    residue_id_alternative : list of str
        Alternative residue ID(s) in case input ID is not available.
    residue_ix : str
        Residue index.
    color : str
        Residue color (matplotlib name).
    _pocket_name : str
        Pocket name.
    _subpocket_name : str
        Subpocket name.
    """

    def __init__(
        self,
        center,
        residue_id,
        residue_id_alternative,
        residue_ix,
        color,
        subpocket_name,
        pocket_name,
    ):
        """
        Initialize an AnchorResidue object.

        Parameters
        ----------
        center : list or numpy.array
            Coordinates (x, y, z) of the residue center.
        residue_id : int
            Residue ID.
        residue_id_alternative : list of int or None
            Alternative residue ID(s) in case input ID is not available.
        residue_ix : int
            Residue index.
        color : str
            Residue color (matplotlib name).
        subpocket_name : str or None
            Subpocket name.
        pocket_name : str or None
            Pocket name.
        """

        if center is not None:
            if len(center) != 3:
                raise ValueError(f"Length of center vector must be 3 but is {len(center)}.")
            center = np.array(center)
        self.center = center
        self.residue_id = residue_id
        self.residue_id_alternative = residue_id_alternative
        self.residue_ix = residue_ix
        self.color = color
        self._pocket_name = pocket_name
        self._subpocket_name = subpocket_name

        if self.residue_id_alternative:
            _logger.info(
                f"Pocket {self._pocket_name} / subpocket {self._subpocket_name}: "
                f"Anchor residue {self.residue_id} (PDB ID) / {self.residue_ix} (index) "
                f"is missing, use {self.residue_id_alternative} instead."
            )
        if self.center is None:
            _logger.info(
                f"Pocket {self._pocket_name} / subpocket {self._subpocket_name}: "
                f"Anchor residue {self.residue_id} (PDB ID) / {self.residue_ix} (index) "
                f"and its neighbors are missing - no anchor residue could be set."
            )

    @property
    def anchor_residue(self):
        """
        Anchor residue attributes.

        Returns
        -------
        pd.Series
            Anchor residue attributes.
        """
        return pd.Series(
            {
                "anchor_residue.color": self.color,
                "anchor_residue.id": self.residue_id,
                "anchor_residue.id_alternative": self.residue_id_alternative,
                "anchor_residue.ix": self.residue_ix,
                "anchor_residue.center": self.center,
            }
        )
