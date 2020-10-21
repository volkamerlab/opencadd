"""
opencadd.structure.pocket.utils

Defines utility functions.
"""


def _format_residue_ids_and_labels(residue_ids, residue_labels):
    """
    Handle input residue IDs and labels: Must be of same length, cast values to string.

    Parameters
    ----------
    residue_ids : list of str
        Pocket residue IDs.
    residue_labels : list of str
        Pocket residue labels.

    Returns
    -------
    tuple of list of str
        Residue IDs, residue labels.
    """

    # If no residue labels are given, create list of None
    # (list length = number of anchor residue IDs)
    if not residue_labels:
        residue_labels = [None] * len(residue_ids)

    # Check if residue IDs and labels are of same length, if not raise error
    if len(residue_ids) != len(residue_labels):
        raise ValueError(f"Number of residue IDs and labels must be of same length.")

    # Cast residue IDs and labels to strings (except None)
    residue_ids = [str(residue_id) if residue_id else residue_id for residue_id in residue_ids]
    residue_labels = [
        str(residue_label) if residue_label else residue_label for residue_label in residue_labels
    ]

    return residue_ids, residue_labels
