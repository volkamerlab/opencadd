"""
This module handles different tasks concerning hydrogen.
"""
from rdkit import Chem

__all__ = ["remove_hydrogens"]


def remove_hydrogens(mol):
    """Removes any hydrogens from the graph of a molecule.

    This is a wrapper around rdkit.Chem.rdmolops.RemoveHs.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The molecule to be modified.

    Returns
    -------
    mol: rdkit.Chem.Mol
        A new molecule with the hydrogens removed.
    """
    mol = Chem.RemoveHs(mol)
    return mol
