"""
opencadd.databases.klifs.remote.coordinates
Utility functions to work with KLIFS data (remote)

Base functions to work with coordinates.
"""

from pathlib import Path

from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem

from ..klifs_client import KLIFS_CLIENT
from .structures import structures_from_structure_ids
from ..utils import COORDINATES_ENTITIES, COORDINATES_INPUT_FORMATS, COORDINATES_OUTPUT_FORMATS
from ..utils import check_entity_format, file_path


def fetch(structure_id, entity='complex', input_format='mol2', output_format='biopandas', compute2d=True):
    """
    Fetch structural data from KLIFS database in different output formats.

    Parameters
    ----------
    structure_id : str
        KLIFS structure ID.
    entity : str
        Structural entity: complex (default), ligand, pocket, or protein.
    input_format : str
        Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).
    output_format : str
        Output format: text, biopandas (default), or rdkit (only for entity=ligand).
    compute2d : bool
        For entity=ligand only. Compute 2D coordinates (default) or keep 3D coordinates.
    """

    check_entity_format(entity, input_format, output_format)

    # Fetch text from KLIFS
    text = _fetch_text(structure_id, entity, input_format)
    if not text:  # TODO Ask Albert why no remote water
        raise ValueError(f'Entity {entity} is not available remotely but we could ask Albert to add this.')

    # Return different output formats
    if output_format == 'text':
        return text

    elif output_format == 'rdkit':
        return _mol2_text_to_rdkit_mol(text, compute2d)

    elif output_format == 'biopandas':
        if input_format == 'mol2':
            return _mol2_text_to_dataframe(text)
        elif input_format == 'pdb':
            return _pdb_text_to_dataframe(text)


def save(structure_id, output_path, entity='complex', input_format='mol2', in_dir=False):
    """
    Save structural data to file.

    Parameters
    ----------
    structure_id : str
        KLIFS structure ID.
    entity : str
        Structural entity: complex (default), ligand, pocket, or protein.
    input_format : str
        Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).
    in_dir : bool
        Save file in KLIFS directory structure (default: False).
    """

    check_entity_format(entity, input_format)
    output_path = Path(output_path)

    # Get metadata
    metadata = structures_from_structure_ids(structure_id).iloc[0]

    # Set up output path (metadata in the form of directory structure or file name)
    output_path = file_path(
        output_path, 
        metadata.species.upper(), 
        metadata.kinase, 
        metadata.pdb, 
        metadata.alt, 
        metadata.chain, 
        entity, 
        input_format, 
        in_dir
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get text
    text = fetch(structure_id, entity, input_format, output_format='text')
    if not text:  # TODO Ask Albert why no remote water
        raise ValueError(f'Entity {entity} is not available remotely but we could ask Albert to add this.')

    # Save text to file
    with open(output_path, 'w') as f:
        f.write(text)


def _fetch_text(structure_id, entity='complex', input_format='mol2'):
    """
    Get structural data content from KLIFS database as string (text).

    Parameters
    ----------
    structure_id : str
        KLIFS structure ID.
    entity : str
        Structural entity: complex (default), ligand, pocket, or protein.
    input_format : str
        Input file format (fetched from KLIFS): mol2 (default) or pdb (only for entity=complex).

    Returns
    -------
    str
        Structural data.
    """
    
    if entity == 'complex' and input_format == 'mol2':
        return KLIFS_CLIENT.Structures.get_structure_get_complex(
            structure_ID=structure_id
        ).response().result
    elif entity == 'complex' and input_format == 'pdb':
        return KLIFS_CLIENT.Structures.get_structure_get_pdb_complex(
            structure_ID=structure_id
        ).response().result
    elif entity == 'ligand' and input_format == 'mol2':
        return KLIFS_CLIENT.Structures.get_structure_get_ligand(
            structure_ID=structure_id
        ).response().result
    elif entity == 'pocket' and input_format == 'mol2':
        return KLIFS_CLIENT.Structures.get_structure_get_pocket(
            structure_ID=structure_id
        ).response().result
    elif entity == 'protein' and input_format == 'mol2':
        return KLIFS_CLIENT.Structures.get_structure_get_protein(
            structure_ID=structure_id
        ).response().result


def _mol2_text_to_dataframe(mol2_text):
    """
    Get structural data from mol2 text.

    Parameters
    ----------
    mol2_text : str
       Mol2 file content from KLIFS database.

    Returns
    -------
    pandas.DataFrame
        Structural data.
    """

    pmol = PandasMol2()

    try:
        mol2_df = pmol._construct_df(
            mol2_text.splitlines(True),
            col_names=[
                'atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge', 'backbone'
            ],
            col_types=[
                int, str, float, float, float, str, int, str, float, str
            ]
        )
    except ValueError:
        mol2_df = pmol._construct_df(
            mol2_text.splitlines(True),
            col_names=[
                'atom_id', 'atom_name', 'x', 'y', 'z', 'atom_type', 'subst_id', 'subst_name', 'charge'
            ],
            col_types=[
                int, str, float, float, float, str, int, str, float
            ]
        )

    return mol2_df


def _mol2_text_to_rdkit_mol(mol2_text, compute2d=True):
    """
    Get structural data from mol2 text.

    Parameters
    ----------
    mol2_text : str
       Mol2 file content from KLIFS database.
    compute2d : bool
        Compute 2D coordinates for ligand (default).

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        Molecule.
    """

    mol = Chem.MolFromMol2Block(mol2_text)

    if compute2d:
        AllChem.Compute2DCoords(mol)

    return mol


def _pdb_text_to_dataframe(pdb_text):
    """
    Get structural data from pdb text.

    Parameters
    ----------
    pdb_text : str
       Pdb file content from KLIFS database.

    Returns
    -------
    dict of pandas.DataFrame
        Structural data
    """

    ppdb = PandasPdb()

    pdb_dict = ppdb._construct_df(
        pdb_text.splitlines(True)
    )

    print(f'Structural data keys: {pdb_dict.keys()}')

    return pdb_dict