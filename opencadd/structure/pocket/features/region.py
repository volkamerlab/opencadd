"""
opencadd.structure.pocket.region

Defines pocket regions.
"""

import pandas as pd


class Region:
    """
    Class defining a region.

    Attributes
    ----------
    name : str
        Region name.
    residue_ids : list of int
        List of residue IDs defining the region.
    residue_ixs : list of int or None
        List of residue indices.
    color : str
        Region color (matplotlib name).
    """

    def __init__(self, name, residue_ids, residue_ixs=None, color="blue"):

        self.name = name
        self.color = color
        self.residue_ids = residue_ids
        self.residue_ixs = residue_ixs

    @property
    def region(self):
        """
        Region attributes.

        Returns
        -------
        pd.DataFrame
            Region attributes (per row one region residue).
        """

        n_residues = len(self.residue_ids)
        return pd.DataFrame(
            {
                "region.name": [self.name] * n_residues,
                "region.color": [self.color] * n_residues,
                "residue.id": self.residue_ids,
                "residue.ix": self.residue_ixs,
            }
        )
