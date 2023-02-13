"""
Query the RESTful Data API and related endpoints of the RCSB Protein Data Bank (PDB) webservers.
"""

# Standard library
from typing import Literal, Union, Optional, Sequence
import gzip
# 3rd party
import numpy as np
# Self
from opencadd._http_request import response_http_request


# General API endpoints
_ROOT_DATA: str = f"https://data.rcsb.org/rest/v1"
_ROOT_FILE: str = f"https://files.rcsb.org"
_ROOT_SEARCH: str = f"https://search.rcsb.org/rcsbsearch/v2/query"
_END_DATA: str = f"{_ROOT_DATA}/core"
_END_SCHEMA: str = f"{_ROOT_DATA}/schema"
_END_HOLDINGS = f"{_ROOT_DATA}/holdings"


def schema(
        schema_type: Literal[
            'entry',
            'polymer_entity',
            'branched_entity',
            'nonpolymer_entity',
            'polymer_entity_instance',
            'branched_entity_instance',
            'nonpolymer_entity_instance',
            'assembly',
            'chem_comp'
        ]
) -> dict:
    """
    Get the data schema for a data type.

    Parameters
    ----------
    schema_type : {
                    'entry', 'polymer_entity', 'branched_entity', 'nonpolymer_entity',
                    'polymer_entity_instance', 'branched_entity_instance',
                    'nonpolymer_entity_instance', 'assembly', 'chem_comp'
                  }
        Type of the data schema.

    Returns
    -------
    dict

    References
    ----------
    * https://data.rcsb.org/index.html#data-schema
    """
    return _data_query(url=f"{_END_SCHEMA}/{schema_type}")


def entry(pdb_id: str) -> dict:
    """
    Description of a PDB entry at the top level of PDB structural hierarchical data organization.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/entry/{pdb_id}")


def entry_pubmed(pdb_id: str) -> dict:
    """
    Description of a PDB entry's primary citation, annotated by PubMed.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/pubmed/{pdb_id}")


def assembly(pdb_id: str, assembly_id: Union[int, str]) -> dict:
    """
    Description of a structural assembly (quaternary structure) in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    assembly_id : int | str
        Assembly ID of the biological assembly candidate in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/assembly/{pdb_id}/{assembly_id}")


def entity_branched(pdb_id: str, entity_id: Union[int, str]) -> dict:
    """
    Description of a branched entity (molecule) in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    entity_id : int | str
        Entity ID of the branched entity in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/branched_entity/{pdb_id}/{entity_id}")


def entity_nonpolymer(pdb_id: str, entity_id: Union[int, str]) -> dict:
    """
    Description of a non-polymer entity (molecule) in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    entity_id : int | str
        Entity ID of the non-polymer entity in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/nonpolymer_entity/{pdb_id}/{entity_id}")


def entity_polymer(pdb_id: str, entity_id: Union[int, str]) -> dict:
    """
    Description of a polymer entity (molecule) in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    entity_id : int | str
        Entity ID of the polymer entity in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/polymer_entity/{pdb_id}/{entity_id}")


def entity_polymer_uniprot(pdb_id: str, entity_id: Union[int, str]) -> dict:
    """
    UniProt annotations for a macromolecular polymer entity (molecule) in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    entity_id : int | str
        Entity ID of the polymer entity in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/uniprot/{pdb_id}/{entity_id}")


def instance_branched(pdb_id: str, asym_id: Union[int, str]) -> dict:
    """
    Description of an instance (chain) of a branched entity in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    asym_id : int | str
        Instance (chain) ID of the branched entity instance in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/branched_entity_instance/{pdb_id}/{asym_id}")


def instance_nonpolymer(pdb_id: str, asym_id: Union[int, str]) -> dict:
    """
    Description of an instance (chain) of a non-polymer entity in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    asym_id : int | str
        Instance (chain) ID of the non-polymer entity instance in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/nonpolymer_entity_instance/{pdb_id}/{asym_id}")


def instance_polymer(pdb_id: str, asym_id: Union[int, str]) -> dict:
    """
    Description of an instance (chain) of a polymer entity in a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    asym_id : int | str
        Instance (chain) ID of the polymer entity instance in the PDB entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/polymer_entity_instance/{pdb_id}/{asym_id}")


def chemical_component(ccd_id: str) -> dict:
    """
    Description of a chemical component (i.e. ligand, small molecule, monomer) in a PDB entry.

    Parameters
    ----------
    ccd_id : str
        CHEM COMP ID of the chemical component. For protein polymer entities, this is the
        three-letter code for the amino acid. For nucleic acid polymer entities, this is the
        one-letter code for the base.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/chemcomp/{ccd_id}")


def chemical_drugbank(ccd_id: str) -> dict:
    """
    Description of a chemical component (i.e. ligand, small molecule, monomer) in a PDB entry,
    annotated by DrugBank.

    Parameters
    ----------
    ccd_id : str
        CHEM COMP ID of the chemical component. For protein polymer entities, this is the
        three-letter code for the amino acid. For nucleic acid polymer entities, this is the
        one-letter code for the base.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/drugbank/{ccd_id}")


def group_entry(group_id: str) -> dict:
    """
    PDB cluster data for entries, based upon a given aggregation method.

    Parameters
    ----------
    group_id : str
        Group ID, e.g. 'Q3Y9I6'.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/entry_groups/{group_id}")


def group_provenance(provenance_id: str) -> dict:
    """
    Aggregation method used to create groups.

    Parameters
    ----------
    provenance_id : str
        Group provenance ID, e.g. 'provenance_sequence_identity'.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/group_provenance/{provenance_id}")


def group_entity(group_id: str) -> dict:
    """
    PDB cluster data for polymer entities, based upon a given aggregation method.

    Parameters
    ----------
    group_id : str
        Group ID, e.g. 'Q3Y9I6'.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/polymer_entity_groups/{group_id}")


def interface(
        pdb_id: str,
        assembly_id: Union[int, str] = 1,
        interface_id: Union[int, str] = 1,
) -> dict:
    """
    Description of a pairwise polymeric interface in an assembly of a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    assembly_id : int | str, optional, default: 1
        Assembly ID of the biological assembly candidate in the PDB entry.
    interface_id : int | str, optional, default: 1
        Interface ID of the pairwise polymeric interface.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_DATA}/interface/{pdb_id}/{assembly_id}/{interface_id}")


def holdings(
        status: Literal["current", "unreleased", "removed"] = "current"
) -> np.ndarray:
    """
    Get all PDB ID holdings data for a specific entry status.

    Parameters
    ----------
    status : Literal["current", "unreleased", "removed"], optional, default: "current"
        Status of PDB entries to retrieve; either all current, all unreleased, or all removed PDB
        IDs.

    Returns
    -------
    numpy.ndarray[ndim=1, dtype="<U4"]
        A one-dimensional array of PDB IDs.
    """
    if status not in ("current", "unreleased", "removed"):
        raise ValueError(f"{status} is not a valid argument for `data`.")
    return np.array(_data_query(url=f"{_END_HOLDINGS}/{status}/entry_ids"))


def holdings_without_pdb_file() -> np.ndarray:
    """
    PDB IDs of all entries without a corresponding legacy PDB-format file.

    Returns
    -------
    numpy.ndarray[ndim: 1, dtype: <U4]
        Array of PDB IDs.

    Notes
    -----
    * Following entries don't have a corresponding PDB-format file:
      * Entries with multiple-character chain IDs.
      * Entries with more than 62 chains.
      * Entries with 100,000 or more atoms.
      * Entries with a complex beta sheet topology.
    * Number of these entries will continue to grow as new large structures are deposited and released.
    * These entries can also be found using Advanced Search (Deposition > Compatible with PDB Format > equals > N)

    References
    ----------
    * `RCSB Documentation: Structures Without Legacy PDB Format Files <https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files>`_
    """
    pdb_ids = response_http_request(
        url="https://files.wwpdb.org/pub/pdb/compatible/pdb_bundle/pdb_bundle_index.txt",
        response_type="str"
    )
    return np.array(pdb_ids.upper().splitlines())


def holdings_status(pdb_id: Union[str, Sequence[str]]) -> dict:
    """
    Status and status code of a PDB entry.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_HOLDINGS}/status/{pdb_id}")


def holdings_removed(pdb_id: str) -> dict:
    """
    Description of an entry that was removed from the PDB repository.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_HOLDINGS}/removed/{pdb_id}")


def holdings_unreleased(pdb_id: str) -> dict:
    """
    Description of an entry that is being processed or on hold waiting for release.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.

    Returns
    -------
    dict
    """
    return _data_query(url=f"{_END_HOLDINGS}/unreleased/{pdb_id}")


def _data_query(url: str) -> dict:
    """
    Send data query and get results in dict format.

    Parameters
    ----------
    url : str
        Full URL of the API query.

    Returns
    -------
    dict
    """
    return response_http_request(url=url, response_type="json")
