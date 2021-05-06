"""
This module removes salt fragments.
"""
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
import csv
import logging
from rdkit import RDLogger
from ...utils import data_path

__all__ = ["remove_salts"]


def _extract_row_from_csv(row_num, filename="salts.tsv", delimiter="\t"):
    """Extracts a row from a  CSV File"""
    with open(data_path(filename)) as csvfile:
        input_csv = csv.reader(csvfile, delimiter=delimiter)
        new_array = []
        for row in input_csv:
            new_read = row[row_num]
            new_array.append(new_read)
    return new_array


def remove_salts(mol, dictionary=True, *args, **kwargs):
    """Removes salts from a molecule.

    This function removes detected salts following a salts dictionary by
    default.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        The molecule to be modified.
    dictionary: bool, optional
        True (default): Activates the use of the salt dictionary.
        False: Uses the standard StripMol functionality, provided by
        rdkit.Chem.SaltRemover.
    defnData: list of str, optional
        If the dictionary is set to False, a custom dictionary can be
        set up. If not rdkit default values from
        '/scratch/RDKit_git/Data/Salts.txt' are used.

    Returns
    -------
    mol: rdkit.Chem.Mol
        A new molecule with salts removed.

    Notes
    -----
    The Salts Dictionary
        The dictionary used is a derived version from the ChEMBL salt
        dictionary, created for the standardiser application by Francis 
        Atkinson. The salts are stored as list of (neutral) SMILES.
    """
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.ERROR)
    i = 0

    if dictionary == True:
        salts = _extract_row_from_csv(0)
        salt_names = _extract_row_from_csv(1)
        list_len = len(salts)

        while i < list_len:
            salt = salts[i]
            salt_name = salt_names[i]
            test = Chem.MolToSmiles(mol)
            i += 1
            remover = SaltRemover(defnData=salt)
            stripped_mol = remover.StripMol(mol)
            if stripped_mol.GetNumAtoms() == 0:
                print(test)
                break
            test_smiles = Chem.MolToSmiles(stripped_mol)
            if test_smiles != test:
                logging.debug("Following salt was stripped: %s", salt_name)
                mol = stripped_mol
    else:
        mol = SaltRemover(*args, **kwargs).StripMol(mol)

    return mol
