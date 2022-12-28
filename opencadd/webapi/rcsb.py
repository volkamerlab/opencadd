"""
Implementation of the RESTful API of Protein Data Bank (PDB) at https://rcsb.org

References
----------
Programmatic access to RCSB:
* https://www.rcsb.org/docs/programmatic-access
Data API:
* https://data.rcsb.org/
RESTful API documentation:
* https://data.rcsb.org/redoc/
Data attributes in responses:
* https://data.rcsb.org/data-attributes.html
File services:
* https://www.rcsb.org/docs/programmatic-access/file-download-services
* https://www.wwpdb.org/ftp/pdb-ftp-sites
"""

# Standard library
from typing import Literal, Union, Optional, Sequence
import gzip
# 3rd party
import numpy as np
# Self
from opencadd.webapi.http_request import response_http_request


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
    schema_type : Literal[
                    'entry', 'polymer_entity', 'branched_entity', 'nonpolymer_entity',
                    'polymer_entity_instance', 'branched_entity_instance',
                    'nonpolymer_entity_instance', 'assembly', 'chem_comp'
                  ]
        Type of the data schema.

    Returns
    -------
    dict

    References
    ----------
    * https://data.rcsb.org/index.html#data-schema
    """
    return _data_query(url=f"{_END_SCHEMA}/{schema_type}")


def data_entry(pdb_id: str) -> dict:
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


def data_entry_pubmed(pdb_id: str) -> dict:
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


def data_assembly(pdb_id: str, assembly_id: Union[int, str]) -> dict:
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


def data_entity_branched(pdb_id: str, entity_id: Union[int, str]) -> dict:
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


def data_entity_nonpolymer(pdb_id: str, entity_id: Union[int, str]) -> dict:
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


def data_entity_polymer(pdb_id: str, entity_id: Union[int, str]) -> dict:
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


def data_entity_polymer_uniprot(pdb_id: str, entity_id: Union[int, str]) -> dict:
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


def data_instance_branched(pdb_id: str, asym_id: Union[int, str]) -> dict:
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


def data_instance_nonpolymer(pdb_id: str, asym_id: Union[int, str]) -> dict:
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


def data_instance_polymer(pdb_id: str, asym_id: Union[int, str]) -> dict:
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


def data_chemical_component(ccd_id: str) -> dict:
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


def data_chemical_drugbank(ccd_id: str) -> dict:
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


def data_group_entry(group_id: str) -> dict:
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


def data_group_provenance(provenance_id: str) -> dict:
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


def data_group_entity(group_id: str) -> dict:
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


def data_interface(
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


def data_holdings(
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


def data_holdings_status(pdb_id: Union[str, Sequence[str]]) -> dict:
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


def data_holdings_removed(pdb_id: str) -> dict:
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


def data_holdings_unreleased(pdb_id: str) -> dict:
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


def file_pdb_entry(
        pdb_id: str,
        file_format: Literal["cif", "pdb", "xml", "bcif"] = "cif",
        biological_assembly_id: Optional[int] = None,
) -> bytes:
    """
    Download a PDB entry file (or a specific biological assembly in an entry file)
    in one of available formats.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    file_format : Literal["cif", "pdb", "xml", "bcif"], optional, default: "cif"
        Format of the file.
    biological_assembly_id : int, optional, default: None
        Assembly ID inside the entry. If provided, only the given assembly will be downloaded,
        otherwise (i.e. when `None`), the whole entry is downloaded.

    Returns
    -------
    bytes
        Content of the downloaded file in bytes.
    """
    # Check whether PDB ID is valid:
    if not isinstance(pdb_id, str):
        raise ValueError("`pdb_id` must be of type string.")
    if len(pdb_id) != 4:
        raise ValueError("`pdb_id` must have 4 characters.")
    if not pdb_id[0].isnumeric() or pdb_id[0] == "0":
        raise ValueError("First character of `pdb_id` must be a non-zero digit.")
    # Build HTTP request URL
    url_prefix = f"{_ROOT_FILE}/download/{pdb_id}"
    if file_format not in ("cif", "pdb", "xml", "bcif"):
        raise ValueError(f"File format {file_format} not recognized.")
    if biological_assembly_id is None:
        url = f"{url_prefix}.{file_format}.gz"
    elif file_format == "cif":
        url = f"{url_prefix}-assembly{biological_assembly_id}.cif.gz"
    elif file_format == "pdb":
        url = f"{url_prefix}.pdb{biological_assembly_id}"
    else:
        raise ValueError("Biological assemblies can only be downloaded in CIF and PDB formats.")
    # Download file, decompress and return as bytes
    return gzip.decompress(response_http_request(url=url, response_type="bytes"))


def file_small_molecule(
        ligand_id: str,
        file_type: Literal["model_coords", "ideal_coords", "def"] = "model_coords",
        file_format: Literal["sdf", "mol2", "cif"] = "sdf",
) -> bytes:
    """
    Download a small molecule file in one of available formats.

    Parameters
    ----------
    ligand_id : str
        ID of the ligand.
    file_type: Literal["model_coords", "ideal_coords", "def"], optional, default: "model_coords"
        Type of the file: model coordinates, ideal coordinates, or definition.
        Notice that coordinates are only available in MOL2 and SDF file formats,
        while definition is only available in CIF. Thus, if the argument of `file_type` is `def`,
        then `file_format` argument is ignored and a CIF file is returned. If the argument of
        `file_type` is `model_coords` or `ideal_coords` and `file_type` is `cif`,
        a `ValueError` is raised.
    file_format : Literal["sdf", "mol2", "cif"], optional, default: "sdf"
        Format of the file.

    Returns
    -------
    bytes
        Content of the downloaded file in bytes.
    """
    if file_format not in ("cif", "mol2", "sdf"):
        raise ValueError(f"File format {file_format} not recognized.")
    if file_type not in ("model_coords", "ideal_coords", "def"):
        raise ValueError(f"File type {file_type} not recognized.")

    _URL_PREFIX_MOL = f"{_ROOT_FILE}/ligands/download/"
    if file_type == "def":
        url = f"{_URL_PREFIX_MOL}{ligand_id}.cif"
    elif file_type == "model_coords":
        if file_format in ("mol2", "sdf"):
            url = f"{_URL_PREFIX_MOL}{ligand_id}_model.{file_format}"
        else:
            raise ValueError(f"File format {file_format} is not accepted for type `model_coords`")
    else:
        if file_format in ("mol2", "sdf"):
            url = f"{_URL_PREFIX_MOL}{ligand_id}_ideal.{file_format}"
        else:
            raise ValueError(f"File format {file_format} is not accepted for type `model_coords`")
    return response_http_request(url=url, response_type="bytes")


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

