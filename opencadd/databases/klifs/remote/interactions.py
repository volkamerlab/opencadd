"""
opencadd.databases.klifs.remote.interactions
Utility functions to work with KLIFS data (remote)

Interaction details.
"""

from ..utils import  _abc_idlist_to_dataframe
from ..klifs_client import KLIFS_CLIENT


def interaction_types():
    """
    Get KLIFS interaction types.

    Returns
    -------
    pandas.DataFrame
        KLIFS interaction types.
    """

    result = KLIFS_CLIENT.Interactions.get_interactions_get_types().response().result

    return _abc_idlist_to_dataframe(result)


def interaction_fingerprint_from_structure_ids(structure_ids):
    """
    Get interaction fingerprint(s) from KLIFS structure ID(s).

    Parameters
    ----------
    structure_ids : int or list of int
        KLIFS structure ID(s).

    Returns
    -------
    pandas.DataFrame
        KLIFS interaction fingerprint(s) (KLIFS structure ID and interaction fingerprint).
    """

    if isinstance(structure_ids, int):
        structure_ids = [structure_ids]

    result = KLIFS_CLIENT.Interactions.get_interactions_get_IFP(
        structure_ID=structure_ids
    ).response().result

    return _abc_idlist_to_dataframe(result)


def klifs_pocket_numbering_from_structure_id(structure_id):
    """
    Get KLIFS pocket numbering (PDB vs. KLIFS numbering).

    Parameters
    ----------
    structure_id : int
        KLIFS structure ID.

    Returns
    -------
    pandas.DataFrame
        KLIFS pocket numbering (KLIFS numbering, PDB numbering and KLIFS position
        (<structural element>.<residue index>)).
    """

    result = KLIFS_CLIENT.Interactions.get_interactions_match_residues(
        structure_ID=structure_id
    ).response().result

    return _abc_idlist_to_dataframe(result)
