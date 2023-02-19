"""
Download various data files from the Protein Data Bank (PDB) webservers.
"""

# Standard library
from typing import Literal, Union, Optional
import gzip
from pathlib import Path
# Self
from opencadd._http_request import response_http_request
from opencadd._typing import PathLike
from opencadd import _helpers_sysio

# General API endpoints
_ROOT_FILE: str = f"https://files.rcsb.org"


def entry(
        pdb_id: str,
        file_format: Literal["cif", "pdb", "xml", "bcif"] = "cif",
        biological_assembly_id: Optional[Union[int, str]] = None,
        output_path: Optional[PathLike] = None,
) -> Union[bytes, Path]:
    """
    Download a PDB entry file in one of available formats.

    Parameters
    ----------
    pdb_id : str
        PDB ID of the entry.
    file_format : {'cif', 'pdb', 'xml', 'bcif'}, optional, default: 'cif'
        Format of the entry file to download.
    biological_assembly_id : int, optional, default: None
        Biological assembly ID of an assembly within the entry.
        If not provided (i.e. when set `None`; default), the asymmetric unit will be downloaded,
        otherwise the file containing the coordinates of the given assembly.
        Notice that many records are only available in the PDB file of the asymmetric unit.
    output_path : str or pathlib.Path
        Path to a local directory for storing the downloaded file.
        If the directory does not exist, it and all its necessary parent directories will be created.
        The filename will be the PDB ID, and the extension will be the same as the `file_format` argument.
        If not provided (i.e. when set `None`; default),
        the byte contents of the downloaded file will be returned.

    Returns
    -------
    bytes or pathlib.Path
        Either the content of the downloaded file in bytes (when `output_path` is `None`),
        or the full filepath of the stored file (when `output_path` is specified).

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
    # Download file and decompress
    byte_content = gzip.decompress(response_http_request(url=url, response_type="bytes"))
    if output_path is None:
        return byte_content
    return _helpers_sysio.save_to_file(content=byte_content, filename=pdb_id, extension=file_format, path=output_path)


def small_molecule(
        ligand_id: str,
        file_type: Literal["model_coords", "ideal_coords", "def"] = "model_coords",
        file_format: Literal["sdf", "mol2", "cif"] = "sdf",
        output_path: Optional[PathLike] = None,
) -> Union[bytes, Path]:
    """
    Download a small molecule file in one of available formats.

    Parameters
    ----------
    ligand_id : str
        ID of the ligand.
    file_type: {'model_coords', 'ideal_coords', 'def'}, optional, default: 'model_coords'
        Type of the file: model coordinates, ideal coordinates, or definition.
        Notice that coordinates are only available in MOL2 and SDF file formats,
        while definition is only available in CIF. Thus, if the input argument is 'def',
        then the `file_format` argument is ignored and a CIF file is returned.
        If the input argument is `model_coords` or `ideal_coords` and `file_type` is set to `cif`,
        a `ValueError` is raised.
    file_format : {'sdf', 'mol2', 'cif'}, optional, default: 'sdf'
        Format of the file to download.
    output_path : str or pathlib.Path
        Path to a local directory for storing the downloaded file.
        If the directory does not exist, it and all its necessary parent directories will be created.
        The filename will be the PDB ID, and the extension will be the same as the `file_format` argument.
        If not provided (i.e. when set `None`; default),
        the byte contents of the downloaded file will be returned.

    Returns
    -------
    bytes or pathlib.Path
        Either the content of the downloaded file in bytes (when `output_path` is `None`),
        or the full filepath of the stored file (when `output_path` is specified).
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

    byte_content = response_http_request(url=url, response_type="bytes")
    if output_path is None:
        return byte_content
    return _helpers_sysio.save_to_file(
        content=byte_content,
        filename=ligand_id,
        extension=file_format,
        path=output_path
    )

