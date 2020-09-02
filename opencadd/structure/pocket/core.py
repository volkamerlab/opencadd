"""
opencadd.structure.pocket.core

Defines core classes and functions.
"""


class Base:
    """
    Basic methods needed for child classes in this module:
    Pocket, Subpocket, Region, AnchorResidue.
    """

    def from_dataframe(self, dataframe, name, color, residue_pdb_ids, residue_labels):
        """
        Set class attributes from DataFrame.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z"
        name : str
            Name.
        color : str
            Color.
        residue_pdb_ids : list of (int, str)
            List of residue PDB IDs.
        residue_labels : list of (int, str) or None
            List of residue labels. Must be of same length as residue_pdb_ids.
        """
        raise NotImplementedError("Implement in your subclass!")

    def _format_residue_pdb_ids_and_labels(self, residue_pdb_ids, residue_labels):
        """
        Handle input residue PDB IDs and labels: Must be of same length, cast values to string.

        Parameters
        ----------

        Returns
        -------
        tuple of list of str
            Residue PDB IDs, residue labels.
        """

        # If no residue labels are given, create list of None
        # (list length = number of anchor residue PDB IDs)
        if not residue_labels:
            residue_labels = [None] * len(residue_pdb_ids)

        # Check if PDB IDs and labels are of same length, if not raise error
        if len(residue_pdb_ids) != len(residue_labels):
            raise ValueError(f"Number of residue PDB IDs and labels must be of same length.")

        # Cast residue PDB IDs and labels to strings (except None)
        residue_pdb_ids = [
            str(residue_pdb_id) if residue_pdb_id else residue_pdb_id
            for residue_pdb_id in residue_pdb_ids
        ]
        residue_labels = [
            str(residue_label) if residue_label else residue_label
            for residue_label in residue_labels
        ]

        return residue_pdb_ids, residue_labels
