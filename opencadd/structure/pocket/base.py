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
        Note: Only residues with int-castable PDB IDs (here: i.e. are not None).

        Returns
        -------
        pandas.DataFrame
            Residue ID and residue index (columns) for all pocket residues (rows).
        """

        residues = {"residue.id": self._residue_ids, "residue.ix": self._residue_ixs}
        residues = (
            pd.DataFrame(residues)
            .dropna(axis=0, subset=["residue.id"])
            .astype({"residue.id": "int"})
        )
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

        # Cast all residue indices to str (raise error if not int-castable)
        if residue_ixs:
            try:
                residue_ixs = [int(residue_ix) for residue_ix in residue_ixs]
            except ValueError as e:
                raise ValueError(
                    f"{e}. Please use only int-castable residue indices. "
                    f"Your input: {residue_ixs}."
                )
            except TypeError as e:
                raise ValueError(
                    f"{e}. Please use only int-castable residue indices. "
                    f"Your input: {residue_ixs}"
                )

        # If no residue indices are given, create list of None
        # (list length = number of anchor residue IDs)
        if not residue_ixs:
            residue_ixs = [None] * len(residue_ids)

        # Check if residue IDs and indices are of same length, if not raise error
        if len(residue_ids) != len(residue_ixs):
            raise ValueError(f"Number of residue IDs and indices must be of same length.")

        # Assign None to residue PDB IDs that cannot be cast to an integer
        residue_ids_clean = []
        residues_cast_to_none = []
        for residue_id, residue_ix in zip(residue_ids, residue_ixs):
            try:
                residue_id_clean = int(residue_id)
            except (ValueError, TypeError):
                residue_id_clean = None
                residues_cast_to_none.append((residue_id, residue_ix))
            residue_ids_clean.append(residue_id_clean)

        if len(residues_cast_to_none) > 0:
            _logger.info(
                f"Pocket {self.name} ({log_text}): "
                f"The following input residues PDB IDs were assigned to the value None "
                f"because they cannot be cast to an integer "
                f"(residue PDB ID, residue index): {residues_cast_to_none}"
            )

        # Raise error if residue IDs or indices contain integer duplicates
        if len(set(residue_ids_clean)) < len(
            [residue_id for residue_id in residue_ids_clean if residue_id]
        ):
            raise ValueError(
                f"Residue PDB IDs cannot contain integer duplicates. Your input: {residue_ids}"
            )
        if residue_ixs and (
            len(set(residue_ixs)) < len([residue_ix for residue_ix in residue_ixs if residue_ix])
        ):
            raise ValueError(
                f"Residue indices cannot contain integer duplicates. Your input: {residue_ixs}"
            )

        return residue_ids_clean, residue_ixs

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
        # Keep only residues that have a residue index assigned
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
