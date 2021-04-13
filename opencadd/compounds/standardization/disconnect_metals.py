"""
This module disconnects metals from structures.
"""
from rdkit.Chem.MolStandardize import rdMolStandardize

__all__ = ["disconnect_metals"]


def disconnect_metals(mol):
    """Disconnects metal atoms that are defined as covalently bonded to
    non-metals.

    Parameters
    ---------
    mol: rdkit.Chem.Mol
        The molecule to be modified.

    Returns
    -------
    mol: rdkit.Chem.Mol
        A new molecule with metals disconnected.
    """
    return rdMolStandardize.MetalDisconnector().Disconnect(mol)
