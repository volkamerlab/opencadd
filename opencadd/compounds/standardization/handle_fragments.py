"""
This module handles different operations with fragments.
"""
from rdkit.Chem.MolStandardize import rdMolStandardize

__all__ = ["remove_fragments", "choose_largest_fragment"]


def remove_fragments(mol):
    """Filters out fragments.

    A predefined list contains numerous known fragments which can be
    filtered out.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule with various fragments.

    Returns
    -------
    mol: rdkit.Chem.Mol
        Returns a molecule filtered from known fragments.

    Notes
    -----
    The predefined list containing fragments is a list REMOVE_FRAGMENTS
    saved in rdkit/Chem/MolStandardize/fragments.py
    """
    return rdMolStandardize.FragmentRemover().remove(mol)


def choose_largest_fragment(mol):
    """Gets the largest fragment.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule with various fragments in various sizes.

    Returns
    -------
    mol: rdkit.Chem.Mol
        Returns a molecule containing the largest fragment.
    """
    return rdMolStandardize.LargestFragmentChooser().choose(mol)
