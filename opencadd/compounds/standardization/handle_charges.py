"""
This module performs changes on charges.
"""
from rdkit.Chem.MolStandardize import rdMolStandardize

__all__ = ["uncharge", "reionize"]


def uncharge(mol):
    """Attempts to neutralize charges by adding and/or removing
    hydrogens where possible.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The molecule where the charges have to be neutralized.

    Returns
    -------
    mol: rdkit.Chem.Mol
        Returns a neutralized molecule.
    """
    return rdMolStandardize.Uncharger().uncharge(mol)


def reionize(mol):
    """Ensure the strongest acid groups ionize first in partially
    ionized molecules.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The partially ionized molecule.

    Returns
    -------
    mol: rdkit.Chem.Mol
        Returns a molecule, with strongest acid groups ionized first.
    """
    return rdMolStandardize.Reionizer().reionize(mol)
