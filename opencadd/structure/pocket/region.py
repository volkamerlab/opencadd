"""
opencadd.structure.pocket.region

Defines pocket regions.
"""

import pandas as pd

from .core import Base


class Region(Base):
    """
    Class defining a region.

    Attributes
    ----------
    name : str
        Region name.
    color : str
        Region color (matplotlib name).
    residue_pdb_ids : list of (int, str)
        List of residue PDB IDs defining the region.
    residue_labels : list of (int, str)
        List of residue labels.
    """

    def __init__(self):

        self.name = None
        self.color = None
        self.residue_pdb_ids = None
        self.residue_labels = None

    def from_dataframe(self, dataframe, name, color, residue_pdb_ids, residue_labels=None):
        """
        Set region properties.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z"
        name : str
            Region name.
        color : str
            Region color.
        residue_pdb_ids : list of (int, str)
            List of residue PDB IDs defining the region.
        residue_labels : list of (int, str) or None
            List of residue labels. Must be of same length as residue_pdb_ids.
        """

        self.name = name
        self.color = color

        # Format residue PDB IDs and labels
        residue_pdb_ids, residue_labels = self._format_residue_pdb_ids_and_labels(
            residue_pdb_ids, residue_labels
        )

        # Add residue labels to dataframe
        residue_labels_df = pd.DataFrame(
            {"residue.pdb_id": residue_pdb_ids, "residue.label": residue_labels}
        )
        dataframe = dataframe.merge(residue_labels_df, on="residue.pdb_id", how="left")

        # Keep only existing residue PDB IDs
        residues = dataframe[["residue.pdb_id", "residue.label"]].drop_duplicates()
        residues.reset_index(drop=True, inplace=True)
        residues = residues[residues["residue.pdb_id"].isin(residue_pdb_ids)]
        self.residue_pdb_ids = residues["residue.pdb_id"].to_list()
        self.residue_labels = residues["residue.label"].to_list()
