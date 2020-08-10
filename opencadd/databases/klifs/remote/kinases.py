"""
opencadd.databases.klifs.remote.kinases
Utility functions to work with KLIFS data (remote)

Kinase details.
"""

import pandas as pd

from ..utils import _abc_idlist_to_dataframe
from ..klifs_client import KLIFS_CLIENT


def kinase_groups():
    """
    Get all kinase groups.

    Returns
    -------
    list of str
        Kinase group names.
    """

    return KLIFS_CLIENT.Information.get_kinase_groups().response().result


def kinase_families(kinase_group=None):
    """
    Get all kinase families for a kinase group.

    Parameters
    ----------
    kinase_group : None or str
        Kinase group name (default is None, i.e. all kinase groups are selected).

    Returns
    -------
    list of str
        Kinase family names.
    """

    return KLIFS_CLIENT.Information.get_kinase_families(
        kinase_group=kinase_group
    ).response().result


def kinase_names(kinase_group=None, kinase_family=None, species=None):
    """
    Get all kinase names for kinases belonging to a given kinase group, kinase family and/or species (default is None,
    i.e. get all kinase names). If multiple parameters are set, only kinases fullfilling all conditions are returned.

    Parameters
    ----------
    kinase_group : None or str
        Kinase group name (default is None, i.e. all kinase groups are selected).
    kinase_family : None or str
        Kinase family name (default is None, i.e. all kinase families are selected).
    species : None or str
        Species name (default is None, i.e. all species are selected).

    Returns
    -------
    pandas.DataFrame
        Kinase names with details (kinase KLIFS ID, kinase name (abbreviation and full), and species).
    """

    results = KLIFS_CLIENT.Information.get_kinase_names(
        kinase_group=kinase_group,
        kinase_family=kinase_family,
        species=species
    ).response().result

    return _abc_idlist_to_dataframe(results)


def kinases_from_kinase_names(kinase_names, species=None):
    """
    Get all kinases (+details) by kinase name(s).

    Parameters
    ----------
    kinase_names : str or list of str
        Kinase name.
    species : None or str
        Species name (default is None, i.e. all species are selected).

    Returns
    -------
    pandas.DataFrame
        Kinase(s) details (kinase KLIFS IDs, names (abbreviation, HGNC, and full), family, gropu, kinase class, species,
        UniProt ID, IUPHAR ID and pocket sequence).
    """

    if isinstance(kinase_names, str):
        kinase_names = [kinase_names]

    results = []

    for kinase_name in kinase_names:

        try:

            result = KLIFS_CLIENT.Information.get_kinase_ID(
                kinase_name=kinase_name,
                species=species
            ).response().result
            result_df = _abc_idlist_to_dataframe(result)
            results.append(result_df)

        except Exception as e:
            print(f'Kinase {kinase_name}: {e}')

    # Kinase IDs can occur multiple times if the input kinase names describe the same kinase, thus drop duplicates
    kinases = pd.concat(results)
    kinases = kinases.drop_duplicates('kinase_ID').reset_index(drop=True)

    return kinases


def kinases_from_kinase_ids(kinase_ids):
    """
    Get all kinases (+details) by KLIFS kinase ID(s).

    Parameters
    ----------
    kinase_ids : int or list of int
        KLIFS kinase ID(s).

    Returns
    -------
    pandas.DataFrame
        Kinase(s) details (kinase KLIFS IDs, names (abbreviation, HGNC, and full), family, gropu, kinase class, species,
        UniProt ID, IUPHAR ID and pocket sequence).
    """

    if isinstance(kinase_ids, int):
        kinase_ids = [kinase_ids]

    results = []

    for kinase_id in kinase_ids:

        result = KLIFS_CLIENT.Information.get_kinase_information(
            kinase_ID=[kinase_id]
        ).response().result
        result_df = _abc_idlist_to_dataframe(result)
        results.append(result_df)

    return pd.concat(results)
