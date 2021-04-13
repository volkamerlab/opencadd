"""
This module performs tautomer enumeration and canonicalization.
"""
from rdkit.Chem.MolStandardize.standardize import (
    enumerate_tautomers_smiles,
    canonicalize_tautomer_smiles,
)

__all__ = ["enumerate_tautomer", "canonicalize_tautomer"]


def enumerate_tautomer(smiles):
    """Generates all possible tautomers.

    During the enumeration it generates all possible tautomers using a
    series of tranformation rules. It also removes stereochemistry from
    double bonds that are single in at least 1 tautomer.

    Parameters
    ----------
    smiles: str

    Returns
    -------
    possible_tautomers: list of str

    Notes
    -----
    The default list of Tautomer Transforms is saved in the list:
    TAUTOMER_TRANSFORMS in rdkit/Chem/MolStandardize/tautomer.py

    """
    possible_tautomers = enumerate_tautomers_smiles(smiles)
    return possible_tautomers


def canonicalize_tautomer(smiles):
    """Generates a canonicalized tautomer.

    During the canonicalization it also enumerates all possible
    tautomers, but after that it uses a scoring system to determine a
    canonical tautomer.

    Parameters
    ----------
    smiles: str

    Returns
    -------
    smiles: str

    Notes
    -----
    The default list of Tautomer Scores is saved in the list:
    TAUTOMER_SCORES in rdkit/Chem/MolStandardize/tautomer.py
    """
    return canonicalize_tautomer_smiles(smiles)
