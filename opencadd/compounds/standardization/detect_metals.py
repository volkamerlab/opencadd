"""
This module detects all metals.
"""
import logging
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.info')

__all__ = ["detect_metals"]


def _validation_smiles(mol):
    """Utility function that converts a mol to SMILES for later validation.
    """
    validation_smiles = Chem.MolToSmiles(mol)
    return validation_smiles


def detect_metals(mol, *args, **kwargs):
    """Detects metals.

    Generates a SMILES out of the entered mol for validation, performs metal 
    disconnection, turns the changed mol into another SMILES 
    and validates it with the first SMILES created.  

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The molecule which has to be searched for non-organic
        substructures.

    Returns
    -------
    boolean: bool
        Returns if the stucture contains a metal (True),
        or not (False).
    """

    smiles_before = _validation_smiles(mol)
    mol_without_metal = rdMolStandardize.MetalDisconnector().Disconnect(mol)
    smiles_after = _validation_smiles(mol_without_metal)
    if smiles_before == smiles_after:
        return False
    else:
        return True
