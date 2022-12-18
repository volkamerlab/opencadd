"""
Functions for downloading files from the PDB archive.

References
----------
https://www.rcsb.org/docs/programmatic-access/file-download-services
https://www.wwpdb.org/ftp/pdb-ftp-sites
"""

import gzip
from typing import Literal, Optional
from opencadd.api.web.web import response_http_request


# General API endpoints and data
_URL_ROOT = "https://files.rcsb.org/"
_URL_PREFIX_PDB = f"{_URL_ROOT}download/"
_URL_PREFIX_MOL = f"{_URL_ROOT}ligands/download/"


def pdb_entry(
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
    url_prefix = f"{_URL_PREFIX_PDB}{pdb_id}"
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
