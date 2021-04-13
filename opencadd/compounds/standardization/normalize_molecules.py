"""
This module performs Normalization transformations
to correct functional groups and recombine charges.
"""
from rdkit.Chem.MolStandardize import rdMolStandardize

__all__ = ["normalize"]


def normalize(mol):
    """Applies a series of Normalization transforms to correct
    functional groups and recombine charges.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    mol: rdkit.Chem.Mol
        Returns a molecule where various Normalization transforms to
        correct functional groups and recombine charges have been
        performed on.

    Notes
    -----
    The Normalization transformations are saved in the list
    NORMALIZATIONS contained in rdkit/Chem/MolStandardize/normalize.py
    """
    return rdMolStandardize.Normalizer().normalize(mol)
