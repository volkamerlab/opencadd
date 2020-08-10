"""
opencadd.databases.klifs.utils
Utility functions to work with KLIFS data

General utility functions
"""

from pathlib import Path

from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

COORDINATES_ENTITIES = ['complex', 'ligand', 'pocket', 'protein', 'water']
COORDINATES_INPUT_FORMATS = ['mol2', 'pdb']
COORDINATES_OUTPUT_FORMATS = ['text', 'biopandas', 'rdkit']


def file_path(path, species, kinase, pdb, alt, chain, entity='complex', format='mol2', in_dir=False):
    """
    Get file path.

    Parameters
    ----------
    path : pathlib.Path or str
        Path to folder for file destination (if in_dir=False) or KLIFS_downlaod folder destination (if in_dir=True).
    species : str
        Species.
    kinase : str
        Kinase name.
    pdb : str
        PDB ID.
    alt : str
        Alternate model ID.
    chain : str
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

    species = species.upper()
    alt = alt.replace('-', '')
    chain = chain.replace('-', '')
    structure = f"{pdb}{f'_alt{alt}' if bool(alt) else ''}{f'_chain{chain}' if bool(chain) else ''}"

    if in_dir:
        path = path / 'KLIFS_download' / species / kinase / structure / f'{entity}.{format}'
    else:
        path = path / f'{species}_{kinase}_{structure}_{entity}.{format}'

    return path


def check_entity_format(entity, input_format, output_format=None):
    """
    Check if entity and format (and their combinations) are available.

    Parameters
    ----------
    entity : str
        Structural entity: complex, ligand, pocket, protein, or water (only in local module).
    input_format : str
        File input format: mol2 or pdb (only for entity=complex).
    output_format : str or None
        Output format: text (only in remote module), biopandas, or rdkit (only for entity=ligand).
    """

    # Check if parameters are valid
    if entity not in COORDINATES_ENTITIES:
        raise ValueError(f'Invalid entity. Select from {", ".join(COORDINATES_ENTITIES)}.')
    if input_format not in COORDINATES_INPUT_FORMATS:
        raise ValueError(f'Invalid input format. Select from {", ".join(COORDINATES_INPUT_FORMATS)}.')
    if output_format:
        if output_format not in COORDINATES_OUTPUT_FORMATS:
            raise ValueError(f'Invalid output format. Select from {", ".join(COORDINATES_OUTPUT_FORMATS)}.')

    # Check if parameter combination is valid
    if input_format == 'pdb' and entity != 'complex':
        raise ValueError(f'Entity {entity} is only available in mol2 format.')
    if output_format:
        if output_format == 'rdkit' and entity != 'ligand':
            raise ValueError(f'Only entity ligand can be fetched as rdkit molecule.')


def _abc_idlist_to_dataframe(abc_idlist):
    """
    Transform ABC IDList object into DataFrame.

    Parameters
    ----------
    abc_idlist : list of acb.IDList
        List of labeled list objects from abstract base classes module.

    Returns
    -------
    pandas.DataFrame
        Table with list labels as column names.
    """

    result = abc_idlist

    keys = list(result[0])

    results_dict = {key: [] for key in keys}

    for result in abc_idlist:
        for key in keys:
            results_dict[key].append(result[key])

    return pd.DataFrame(results_dict)