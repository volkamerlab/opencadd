"""
This module performs sanitizations on molecules.
"""
from rdkit import Chem

__all__ = ["sanitize_all"]


def sanitize_all(mol, *args, **kwargs):
    """Performs sanitization on molecules.

    Parameters
    ---------
    mol: rdkit.Chem.Mol
        The molecule that has to be sanitized.

    Returns
    -------
    mol: rdkit.Chem.Mol
        A new sanitized molecule.
    """
    return Chem.SanitizeMol(mol, *args, **kwargs)
