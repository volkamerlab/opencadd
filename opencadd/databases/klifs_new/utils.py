"""
utils.py

Defines utility functions.
"""

from contextlib import contextmanager
import logging
from pathlib import Path
import re

from bravado.client import SwaggerClient

_logger = logging.getLogger(__name__)

KLIFS_API_DEFINITIONS = "http://klifs.vu-compmedchem.nl/swagger/swagger.json"
KLIFS_CLIENT = SwaggerClient.from_url(
    KLIFS_API_DEFINITIONS, config={"validate_responses": False}
)


def metadata_to_filepath(  # Rename function metadata_to_filepath TODO
    path_to_klifs_download,
    species,
    kinase_name,
    structure_pdb,
    structure_alternate_model,
    structure_chain,
    entity="complex",
    input_format="mol2",
    in_dir=False,
):
    """
    Get file path.

    Parameters
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder for file destination (if in_dir=False) or KLIFS_download folder destination (if in_dir=True).
    species : str
        Species.
    kinase_name : str
        Kinase name.
    structure_pdb : str
        PDB ID.
    structure_alternate_model : str
        Alternate model ID.
    structure_chain : str
        Chain ID.
    entity : str
        Structural entity: complex (default), ligand, pocket, protein, or water (only in local module).
    input_format : str
        File format: mol2 (default) or pdb (only for entity=complex).
    in_dir : bool
        Use KLIFS directory structure (default: False).

    Returns
    -------
    pathlib.Path
        File path.
    """

    path_to_klifs_download = Path(path_to_klifs_download)
    species = species.upper()
    if not structure_alternate_model:
        structure_alternate_model = "-"
    structure_alternate_model = structure_alternate_model.replace("-", "")
    if not structure_chain:
        structure_chain = "-"
    structure_chain = structure_chain.replace("-", "")
    structure = f"{structure_pdb}{f'_alt{structure_alternate_model}' if bool(structure_alternate_model) else ''}{f'_chain{structure_chain}' if bool(structure_chain) else ''}"

    if in_dir:
        path = (
            path_to_klifs_download
            / species
            / kinase_name
            / structure
            / f"{entity}.{input_format}"
        )
    else:
        path = (
            path_to_klifs_download
            / f"{species}_{kinase_name}_{structure}_{entity}.{input_format}"
        )

    return path


def filepath_to_metadata(filepath):
    """
    Write docs TODO
    """

    # Cast to string
    filepath = str(filepath)

    # Split filepath
    metadata = re.split(r"/|_|\.", filepath)
    # Assign list elements to detail type
    if ("chain" in filepath) and ("alt" in filepath):
        metadata = {
            "species": metadata[-7].capitalize(),
            "kinase_name": metadata[-6],
            "structure_pdb": metadata[-5],
            "structure_alternate_model": metadata[-4][-1],
            "structure_chain": metadata[-3][-1],
            "entity": metadata[-2],
            "input_format": metadata[-1],
        }
    elif ("chain" in filepath) and ("alt" not in filepath):
        metadata = {
            "species": metadata[-6].capitalize(),
            "kinase_name": metadata[-5],
            "structure_pdb": metadata[-4],
            "structure_alternate_model": None,
            "structure_chain": metadata[-3][-1],
            "entity": metadata[-2],
            "input_format": metadata[-1],
        }
    else:
        metadata = {
            "species": metadata[-5].capitalize(),
            "kinase_name": metadata[-4],
            "structure_pdb": metadata[-3],
            "structure_alternate_model": None,
            "structure_chain": None,
            "entity": metadata[-2],
            "input_format": metadata[-1],
        }
    return metadata


@contextmanager
def silence_logging(highest_level=logging.CRITICAL):
    """
    A context manager that will prevent any logging messages
    triggered during the body from being processed.
    
    Parameters
    ----------
    highest_level : 
        The maximum logging level in use.
        This would only need to be changed if a custom level greater than CRITICAL is defined.

    Notes
    -----
    Two kind-of hacks here:
    - Can't get the highest logging level in effect => delegate to the user.
    - Can't get the current module-level override => use an undocumented (but non-private!) 
      interface.

    References
    ----------
    Code taken from: https://gist.github.com/simon-weber/7853144
    """

    previous_level = logging.root.manager.disable

    logging.disable(highest_level)

    try:
        yield
    finally:
        logging.disable(previous_level)
