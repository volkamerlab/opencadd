"""
opencadd.databases.klifs.remote.ligands
Utility functions to work with KLIFS data (remote)

Ligand details.
"""

import pandas as pd

from ..utils import _abc_idlist_to_dataframe
from ..klifs_client import KLIFS_CLIENT


def ligands_from_kinase_ids(kinase_ids):
    """
    Get ligand ID(s) and details by KLIFS kinase ID(s).

    Parameters
    ----------
    kinase_ids : int or list of int
        KLIFS kinase ID(s).

    Returns
    -------
    pandas.DataFrame
        Ligand(s) details (KLIFS kinase IDs, KLIFS ligand IDs, ligand PDB code, ligand name, SMILES, and InChIKey).
    """

    if isinstance(kinase_ids, int):
        kinase_ids = [kinase_ids]

    results = []

    for kinase_id in kinase_ids:
        result = KLIFS_CLIENT.Ligands.get_ligands_list(kinase_ID=[kinase_id]).response().result
        result_df = _abc_idlist_to_dataframe(result)
        result_df.insert(0, "kinase_id", kinase_id, True)
        results.append(result_df)

    return pd.concat(results)


def structures_from_ligand_ids(ligand_ids):
    """
    Get structure ID(s) and details by KLIFS ligand ID(s).

    Parameters
    ----------
    ligand_ids : int or list of int
        KLIFS ligand ID(s).

    Returns
    -------
    pandas.DataFrame
        Structure(s) details (KLIFS ligand IDs, KLIFS structure IDs, KLIFS kinase IDs, kinase names, species, PDB ID,
        alternate model, chain, RMSD 1+2, pocket sequence, resolution, quality score, missing residues, missing atoms,
        orthosteric and allosteric ligand, DFG conformation, aC helix conformation, G-rich distance, G-rich angle,
        and subpockets).
    """

    if isinstance(ligand_ids, int):
        ligand_ids = [ligand_ids]

    results = []

    for ligand_id in ligand_ids:

        result = (
            KLIFS_CLIENT.Ligands.get_ligands_list_structures(ligand_ID=[ligand_id])
            .response()
            .result
        )
        result_df = _abc_idlist_to_dataframe(result)
        result_df.insert(0, "ligand_id", ligand_id, True)
        results.append(result_df)

    return pd.concat(results)
