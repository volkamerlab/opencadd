"""
opencadd.databases.klifs.remote.structures
Utility functions to work with KLIFS data (remote)

Structure details.
"""

from ..utils import _abc_idlist_to_dataframe
from ..klifs_client import KLIFS_CLIENT


def structures_from_structure_ids(structure_ids):
    """
    Get structure details by KLIFS structure ID(s).

    Parameters
    ----------
    structure_ids : int or list of int
        KLIFS structure ID(s).

    Returns
    -------
    pandas.DataFrame
        Structure details (KLIFS structure IDs, kinase names, species, PDB ID,
        alternate model, chain, RMSD 1+2, pocket sequence, resolution, quality score, missing residues, missing atoms,
        orthosteric and allosteric ligand, DFG conformation, aC helix conformation, G-rich distance, G-rich angle,
        and subpockets).
    """

    if isinstance(structure_ids, int):
        structure_ids = [structure_ids]

    result = (
        KLIFS_CLIENT.Structures.get_structure_list(structure_ID=structure_ids).response().result
    )
    result_df = _abc_idlist_to_dataframe(result)

    return result_df


def structures_from_kinase_ids(kinase_ids):
    """
    Get structure details by KLIFS kinase ID(s).

    Parameters
    ----------
    kinase_ids : int or list of int
        KLIFS kinase ID(s).

    Returns
    -------
    pandas.DataFrame
        Structure details (KLIFS structure IDs, kinase names, species, PDB ID,
        alternate model, chain, RMSD 1+2, pocket sequence, resolution, quality score, missing residues, missing atoms,
        orthosteric and allosteric ligand, DFG conformation, aC helix conformation, G-rich distance, G-rich angle,
        and subpockets).
    """

    if isinstance(kinase_ids, int):
        kinase_ids = [kinase_ids]

    result = KLIFS_CLIENT.Structures.get_structures_list(kinase_ID=kinase_ids).response().result
    result_df = _abc_idlist_to_dataframe(result)

    return result_df


def structures_from_pdb_ids(pdb_ids, alt=None, chain=None):
    """
    Get structure details by PDB ID(s), optionally filter by alternate model and chain.

    Parameters
    ----------
    pdb_ids : str or list of str
        PDB ID(s).
    alt : None or str
        Alternate model. Will only have an effect if only one PDB ID is given.
    chain : None or str
        Chain. Will only have an effect if only one PDB ID is given.

    Returns
    -------
    pandas.DataFrame
        Structure details (KLIFS structure IDs, kinase names, species, PDB ID,
        alternate model, chain, RMSD 1+2, pocket sequence, resolution, quality score, missing residues, missing atoms,
        orthosteric and allosteric ligand, DFG conformation, aC helix conformation, G-rich distance, G-rich angle,
        and subpockets).
    """

    # Ignore alternate model and chain input (set to None) in these two cases:
    # - If multiple PDB IDs are given
    # - If alternate model or chain indicate that no alternate model or chain is selected
    if alt in ["", " ", "-", "_"] or isinstance(pdb_ids, list):
        alt = None
    if chain in ["", " ", "-", "_"] or isinstance(pdb_ids, list):
        chain = None

    if isinstance(pdb_ids, str):
        pdb_ids = [pdb_ids]

    result = KLIFS_CLIENT.Structures.get_structures_pdb_list(pdb_codes=pdb_ids).response().result
    result_df = _abc_idlist_to_dataframe(result)

    # If alt and/or chain are given, filter dataset, else return full dataset
    if alt is not None and chain is None:
        return result_df[result_df.alt == alt]
    elif alt is None and chain is not None:
        return result_df[result_df.chain == chain]
    elif alt is not None and chain is not None:
        return result_df[(result_df.alt == alt) & (result_df.chain == chain)]
    else:
        return result_df
