"""
utils.py

Defines utility functions.
"""

import logging
from pathlib import Path

_logger = logging.getLogger(__name__)


def get_file_path(
    path_to_klifs_download,
    species,
    kinase_name,
    structure_pdb,
    structure_alternate_model,
    structure_chain,
    entity="complex",
    format="mol2",
    in_dir=False,
):
    """
    Get file path.

    Parameters
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder for file destination (if in_dir=False) or KLIFS_downlaod folder destination (if in_dir=True).
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
    entitiy : str
        Structural entity: complex (default), ligand, pocket, protein, or water (only in local module).
    format : str
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
    structure_alternate_model = structure_alternate_model.replace("-", "")
    structure_chain = structure_chain.replace("-", "")
    structure = f"{structure_pdb}{f'_alt{structure_alternate_model}' if bool(structure_alternate_model) else ''}{f'_chain{structure_chain}' if bool(structure_chain) else ''}"

    if in_dir:
        path = (
            path_to_klifs_download
            / "KLIFS_download"
            / species
            / kinase_name
            / structure
            / f"{entity}.{format}"
        )
    else:
        path = (
            path_to_klifs_download
            / f"{species}_{kinase_name}_{structure}_{entity}.{format}"
        )

    return path
