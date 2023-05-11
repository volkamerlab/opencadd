"""
This method assigns stereochemistry to molecules.
"""
from rdkit import Chem

__all__ = ["assign_stereochemistry"]


def assign_stereochemistry(mol, *args, **kwargs):
    """Does Stereochemistry assignment following the Cahn–Ingold–Prelog
    priority rules.

    Parameters
    ---------
    mol: rdkit.Chem.Mol
        The molecule that needs a stereochemistry assignment.

    Returns
    -------
    mol: rdkit.Chem.Mol
        A new molecule with a stereochemistry assigned.
    """
    return Chem.AssignStereochemistry(mol, force=True, cleanIt=True, *args, **kwargs)
