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
    """

    def __init__(self, center, residue_id, residue_id_alternative, residue_ix, color):
        """
        Initialize an AnchorResidue object.

        Parameters
        ----------
        center : list or numpy.array
            Coordinates (x, y, z) of the residue center.
        residue_id : str
            Residue ID.
        residue_id_alternative : list of str or None
            Alternative residue ID(s) in case input ID is not available.
        residue_ix : str
            Residue index.
        color : str
            Residue color (matplotlib name).
        """

        if center:
            if len(center) != 3:
                raise ValueError(f"Length of center vector must be 3 but is {len(center)}.")
            center = np.array(center)
        self.center = center
        self.residue_id = str(residue_id)
        self.residue_id_alternative = residue_id_alternative
        self.residue_ix = str(residue_ix)
        self.color = color

        if self.residue_id_alternative:
            _logger.info(
                f"Anchor residue {self.residue_id} is missing, use {self.residue_id_alternative} instead."
            )
        if self.center is None:
            _logger.info(
                f"Anchor residue {self.residue_id} and its neighbors are missing "
                f"- no anchor residue could be set."
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
