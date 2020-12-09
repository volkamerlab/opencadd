"""
opencadd.structure.pocket.utils

Defines utility functions.
"""

import logging

_logger = logging.getLogger(__name__)


def _format_residue_ids_and_ixs(residue_ids, residue_ixs):
    """
    Handle input residue IDs and indices: Must be of same length, cast values to string.

    Parameters
    ----------
    residue_ids : list of str
        Pocket residue IDs.
    residue_ixs : list of str
        Pocket residue indices.

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

    # Cast residue IDs and indices to strings (except None) TODO
    residue_ids = [str(residue_id) if residue_id else residue_id for residue_id in residue_ids]
    residue_ixs = [str(residue_ix) if residue_ix else residue_ix for residue_ix in residue_ixs]

    # Keep only residues that have IDs/indices that can be cast to an integer
    residue_ids_kept = []
    residue_ixs_kept = []
    residues_dropped = []
    for residue_id, residue_ix in zip(residue_ids, residue_ixs):
        if residue_id:
            try:
                int(residue_id)
                if residue_ix:
                    int(residue_ix)
                residue_ids_kept.append(residue_id)
                residue_ixs_kept.append(residue_ix)
            except ValueError:
                residues_dropped.append((residue_id, residue_ix))
        else:
            residues_dropped.append((residue_id, residue_ix))

    if len(residue_ids) > len(residue_ids_kept):
        _logger.info(
            f"Residues were dropped because they cannot be cast to an integer "
            f"(residue PDB ID, residue index):\n{residues_dropped}"
        )

    return residue_ids_kept, residue_ixs_kept
