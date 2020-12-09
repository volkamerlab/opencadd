"""
opencadd.structure.pocket.region

Defines pocket regions.
"""

import pandas as pd

from .utils import _format_residue_ids_and_ixs


class Region:
    """
    Class defining a region.

    Attributes
    ----------
    name : str
        Region name.
    color : str
        Region color (matplotlib name).
    residue_ids : list of (int, str)
        List of residue IDs defining the region.
    residue_ixs : list of (int, str)
        List of residue indices.
    """

    def __init__(self):

        self.name = None
        self.color = None
        self.residue_ids = None
        self.residue_ixs = None

    def from_dataframe(self, dataframe, name, residue_ids, color="blue", residue_ixs=None):
        """
        Set region properties.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"
        name : str
            Region name.
        residue_ids : list of (int, str)
            List of residue IDs defining the region.
        color : str
            Region color (matplotlib name), blue by default
        residue_ixs : list of (int, str) or None
            List of residue indices. Must be of same length as residue_ids.
        """

        self.name = name
        self.color = color

        # Format residue IDs and indices
        residue_ids, residue_ixs = _format_residue_ids_and_ixs(residue_ids, residue_ixs)

        # Add residue indices to dataframe
        residue_ixs_df = pd.DataFrame({"residue.id": residue_ids, "residue.ix": residue_ixs})
        dataframe = dataframe.merge(residue_ixs_df, on="residue.id", how="left")

        # Keep only existing residue IDs
        residues = dataframe[["residue.id", "residue.ix"]].drop_duplicates()
        residues.reset_index(drop=True, inplace=True)
        residues = residues[residues["residue.id"].isin(residue_ids)]
        self.residue_ids = residues["residue.id"].to_list()
        self.residue_ixs = residues["residue.ix"].to_list()
