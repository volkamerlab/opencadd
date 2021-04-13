"""
This module detects all inorganic substructures.
"""
import logging
from rdkit import Chem

__all__ = ["detect_inorganic"]

inorganic_elements = Chem.MolFromSmarts(
    "[!#1&!#6&!#7&!#8&!#9&!#15&!#16&!#17&!#35&!#53]"
)


def detect_inorganic(mol, *args, **kwargs):
    """Detects all inorganic substructures.

    Has a list of SMARTS which explicitly exludes all organic elemtents
    and searches if there is a substructure match with an non-organic
    element.

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
    To organic structures we count the following:Hydrogen, Carbon,
    Nitrogen, Oxygen, Fluorine, Phosphorus, Sulfur, Chlorine, Bromine,
    Iodine.

    """
    if mol.GetSubstructMatch(inorganic_elements, *args, **kwargs):
        logging.debug("Structure contains a inorganic element")
        return True
    return False
