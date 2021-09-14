"""
opencadd.databases.klifs.utils

Defines utility functions.
"""

from contextlib import contextmanager
import logging
from pathlib import Path
import re

from bravado.client import SwaggerClient

_logger = logging.getLogger(__name__)

KLIFS_API_DEFINITIONS = "https://dev.klifs.net/swagger_v2/swagger.json"
KLIFS_CLIENT = SwaggerClient.from_url(KLIFS_API_DEFINITIONS, config={"validate_responses": False})

PATH_DATA = Path(__file__).parent / ".." / ".." / "data"


def metadata_to_filepath(
    path_to_klifs_download,
    species,
    kinase_name,
    structure_pdb,
    structure_alternate_model,
    structure_chain,
    entity="complex",
    extension="mol2",
    in_dir=False,
):
    """
    Get file path from metadata.

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
    extension : str
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

    # FIXME: The PDB download for ligands in KLIFS is named "klifs_ligand.pdb"
    # instead of "ligand.expo_id". For the time being (until KLIFS maybe streamlines the file name
    # with all the other file names), rename the file here.
    if entity == "ligand" and extension == "pdb":
        entity = "klifs_ligand"

    if in_dir:
        path = path_to_klifs_download / species / kinase_name / structure / f"{entity}.{extension}"
    else:
        path = path_to_klifs_download / f"{species}_{kinase_name}_{structure}_{entity}.{extension}"

    return path


def filepath_to_metadata(filepath):
    """
    Get metadata from file path.

    Parameters
    ----------
    filepath : pathlib.Path
        File path.

    Returns
    -------
    dict
        Metadata in the form of a dictionary with the following keys and values:
        species : str
            Species.
        kinase_name : str
            Kinase name.
        structure_pdb : str
            PDB ID.
        structure_alternate_model : str
            Alternate model ID. None if not existing.
        structure_chain : str
            Chain ID.
        entity : str
            Structural entity: complex, ligand, pocket, protein, or water.
        extension : str
            File format: mol2 or pdb.
    """

    # Cast to string
    filepath = str(filepath)
    # FIXME: The PDB download for ligands in KLIFS is named "klifs_ligand.pdb"
    # instead of "ligand.expo_id". For the time being (until KLIFS maybe streamlines the file name
    # with all the other file names), rename the file here.
    if "klifs_ligand" in filepath:
        filepath = filepath.replace("klifs_ligand", "ligand")

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
            "extension": metadata[-1],
        }
    elif ("chain" in filepath) and ("alt" not in filepath):
        metadata = {
            "species": metadata[-6].capitalize(),
            "kinase_name": metadata[-5],
            "structure_pdb": metadata[-4],
            "structure_alternate_model": None,
            "structure_chain": metadata[-3][-1],
            "entity": metadata[-2],
            "extension": metadata[-1],
        }
    else:
        metadata = {
            "species": metadata[-5].capitalize(),
            "kinase_name": metadata[-4],
            "structure_pdb": metadata[-3],
            "structure_alternate_model": None,
            "structure_chain": None,
            "entity": metadata[-2],
            "extension": metadata[-1],
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
