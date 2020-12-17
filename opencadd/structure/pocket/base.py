"""
opencadd.structure.pocket.base

Defines base pocket.
"""

import logging

import pandas as pd

_logger = logging.getLogger(__name__)


class BasePocket:
    """
    Class defining a base pocket.
    Mainly focusing on pocket residue related functionalities.
    """

    def __init__(self):

        self.name = None
        self._residue_ids = None
        self._residue_ixs = None

    @property
    def residues(self):
        """
        All pocket's residues.

        Returns
        -------
        pandas.DataFrame
            Residue ID and residue index (columns) for all pocket residues (rows).
        """

        residues = {"residue.id": self._residue_ids, "residue.ix": self._residue_ixs}
        residues = pd.DataFrame(residues)
        return residues.reset_index(drop=True)

    @property
    def center(self):
        """
        Pocket center, i.e. the centroid of all input residues' CA atoms.

        Returns
        ----------
        numpy.array
            Pocket center (coordinates).
        """
        raise NotImplementedError("Implement in your subclass!")

    def _format_residue_ids_and_ixs(self, residue_ids, residue_ixs, log_text):
        """
        Handle input residue IDs and indices: Must be of same length, cast values to string.

        Parameters
        ----------
        residue_ids : list of int
            Pocket residue IDs.
            Strings will be case to integers - if not possible residue will be removed!
        residue_ixs : list of (int or None)
            Pocket residue indices.
        log_text : str
            Add information to be logged.

        Returns
        -------
        tuple of list of str
            Residue IDs, residue indices.
        """

        # If no residue indices are given, create list of None
        # (list length = number of anchor residue IDs)
        if not residue_ixs:
            residue_ixs = [None] * len(residue_ids)

        # Check if residue IDs and indices are of same length, if not raise error
        if len(residue_ids) != len(residue_ixs):
            raise ValueError(f"Number of residue IDs and indices must be of same length.")

        # Keep only residues that have IDs/indices that can be cast to an integer
        residue_ids_kept = []
        residue_ixs_kept = []
        residues_dropped = []
        for residue_id, residue_ix in zip(residue_ids, residue_ixs):
            try:
                residue_id = int(residue_id)
                residue_ids_kept.append(residue_id)
                residue_ixs_kept.append(residue_ix)
            except (ValueError, TypeError):
                residues_dropped.append((residue_id, residue_ix))

        if len(residue_ids) > len(residue_ids_kept):
            _logger.info(
                f"Pocket {self.name} ({log_text}): "
                f"The following input residues were dropped because they cannot be cast to an "
                f"integer (residue PDB ID, residue index): {residues_dropped}"
            )

        return residue_ids_kept, residue_ixs_kept

    def _residue_ix2id(self, residue_ix):
        """
        Get residue PDB ID from residue index.

        Parameters
        ----------
        residue_ix : int
            Residue index.

        Returns
        -------
        str
            Residue PDB ID.
        """

        residues = self.residues
        residues = residues[~residues["residue.ix"].isna()]
        try:
            residue_id = residues.set_index("residue.ix").squeeze().loc[residue_ix]
        except KeyError:
            residue_id = None
        return residue_id

    def _residue_id2ix(self, residue_id):
        """
        Get residue index from residue PDB ID.

        Parameters
        ----------
        residue_id : int
            Residue PDB ID.

        Returns
        -------
        str
            Residue index.
        """

        residues = self.residues
        try:
            residue_ix = residues.set_index("residue.id").squeeze().loc[residue_id]
        except KeyError:
            residue_ix = None
        return residue_ix
