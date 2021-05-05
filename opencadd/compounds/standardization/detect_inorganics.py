"""
This module detects all inorganic substructures.
"""
import logging
from rdkit import Chem

__all__ = ["detect_inorganic", "detect_carbon"]

inorganic_elements = Chem.MolFromSmarts(
    "[!#1&!#6&!#7&!#8&!#9&!#15&!#16&!#17&!#34&!#35&!#53]"
)


def detect_inorganic(mol, inorganic_elements=inorganic_elements, *args, **kwargs):
    """Detects all inorganic substructures.

    Has a list of SMARTS which explicitly exludes all elemtents that can occur
    in an organic molecule and searches if there is a substructure match with 
    element not occuring in organic compounds.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The molecule which has to be searched for non-organic
        substructures.

    Returns
    -------
    boolean: bool
        Returns if the stucture contains a non-organic element (True),
        or not (False).

    Notes
    -----
    To organic structures we count the following: Hydrogen, Carbon,
    Nitrogen, Oxygen, Fluorine, Phosphorus, Sulfur, Chlorine, Selenium, Bromine,
    Iodine.

    """
    if mol.GetSubstructMatch(inorganic_elements, *args, **kwargs):
        logging.debug("Structure contains a inorganic element")
        return True
    return False


def detect_carbon(mol, *args, **kwargs):
    """Detects all occurences of Carbon to determine if a compound is organic.

    Checks for the presence of Carbon in the molecule provided. 

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The molecule which has to be searched for Carbon.

    Returns
    -------
    boolean: bool
        Returns if the stucture contains Carbon (True),
        or not (False).
    """
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 6:
            return True
    return False
