"""
opencadd.databases.klifs.local.coordinates
Utility functions to work with KLIFS data (local)

Base functions to work with coordinates.
"""

from pathlib import Path

from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from ..utils import COORDINATES_ENTITIES, COORDINATES_INPUT_FORMATS
from ..utils import check_entity_format


def load(file_path, output_format="biopandas", compute2d=True):
    """
    Load structural data from KLIFS file in different output formats.

    Parameters
    ----------
    file_path : pathlib.Path or str
        Path to KLIFS file.
    entity : str
        Structural entity: complex (default), ligand, pocket, protein, or water.
    compute2d : bool
        For entity=ligand only. Compute 2D coordinates (default) or keep 3D coordinates.
    """

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File does not exist: {file_path}.")

    # Check if parameters are valid
    entity = file_path.stem
    input_format = file_path.suffix[1:]
    check_entity_format(entity, input_format, output_format)

    # Return different output formats
    if output_format == "rdkit":
        return _mol2_file_to_rdkit_mol(str(file_path), compute2d)

    elif output_format == "biopandas":
        if input_format == "mol2":
            return _mol2_file_to_dataframe(file_path).df
        elif input_format == "pdb":
            pass


def _mol2_file_to_dataframe(mol2_file):
    """
    Get structural data from mol2 file.

    Parameters
    ----------
    mol2_file : pathlib.Path or str
       Path to mol2 file.

    Returns
    -------
    pandas.DataFrame
        Structural data.
    """

    mol2_file = Path(mol2_file)

    pmol = PandasMol2()

    try:
        mol2_df = pmol.read_mol2(
            str(mol2_file),
            columns={
                0: ("atom_id", int),
                1: ("atom_name", str),
                2: ("x", float),
                3: ("y", float),
                4: ("z", float),
                5: ("atom_type", str),
                6: ("subst_id", int),
                7: ("subst_name", str),
                8: ("charge", float),
                9: ("backbone", str),
            },
        )

    except ValueError:
        mol2_df = pmol.read_mol2(str(mol2_file))

    return mol2_df


def _mol2_file_to_rdkit_mol(mol2_file, compute2d=True):
    """
    Get structural data from mol2 file.

    Parameters
    ----------
    mol2_file : pathlib.Path or str
       Path to mol2 file.
    compute2d : bool
        Compute 2D coordinates for ligand (default).

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        Molecule.
    """

    mol = Chem.MolFromMol2File(mol2_file)

    if compute2d:
        AllChem.Compute2DCoords(mol)

    return mol
