"""
Implementation of the RESTful API of the Protein Data Bank (PDB) at https://rcsb.org

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
from typing import Literal, Union, Optional
import gzip
# Self
from opencadd._http_request import response_http_request


# General API endpoints
_ROOT_FILE: str = f"https://files.rcsb.org"


def pdb_entry(
        pdb_id: str,
        file_format: Literal["cif", "pdb", "xml", "bcif"] = "cif",
        biological_assembly_id: Optional[Union[int, str]] = None,
) -> bytes:
    """
    Download a PDB entry file (asymmetric unit or a biological assembly)
    in one of available formats.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    file_format : Literal["cif", "pdb", "xml", "bcif"], optional, default: "cif"
        Format of the file.
    biological_assembly_id : int, optional, default: None
        Biological assembly ID. If not provided (i.e. when `None`), the asymmetric unit will be
        downloaded. Notice that many records are only available in the PDB file of the asymmetric unit.

    Returns
    -------
    bytes
        Content of the downloaded file in bytes.

    Notes
    -----
    * Following entries don't have a corresponding PDB-format file:
      * Entries with multiple-character chain IDs.
      * Entries with more than 62 chains.
      * Entries with 100,000 or more atoms.
      * Entries with a complex beta sheet topology.

    References
    ----------
    * `RCSB Documentation: Structures Without Legacy PDB Format Files <https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files>`_
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
        url = f"{url_prefix}.pdb{biological_assembly_id}.gz"
    else:
        raise ValueError("Biological assemblies can only be downloaded in CIF and PDB formats.")
    # Download file, decompress and return as bytes
    return gzip.decompress(response_http_request(url=url, response_type="bytes"))


def small_molecule(
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
    file_type: {'model_coords', 'ideal_coords', 'def'}, optional, default: 'model_coords'
        Type of the file: model coordinates, ideal coordinates, or definition.
        Notice that coordinates are only available in MOL2 and SDF file formats,
        while definition is only available in CIF. Thus, if the argument of `file_type` is `def`,
        then `file_format` argument is ignored and a CIF file is returned. If the argument of
        `file_type` is `model_coords` or `ideal_coords` and `file_type` is `cif`,
        a `ValueError` is raised.
    file_format : {'sdf', 'mol2', 'cif'}, optional, default: 'sdf'
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
