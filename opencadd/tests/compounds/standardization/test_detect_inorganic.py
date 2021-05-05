"""
test for the module `detect_inorganics`
"""
import pytest
import sys
import rdkit

from rdkit import Chem

from opencadd.compounds.standardization import detect_inorganics


def _evaluation_mol_generator(test_smiles=None, test_inchi=None):
    """Creates mol files directly with rdkits functions for evaluation."""
    if test_smiles is not None:
        return Chem.MolFromSmiles(test_smiles)
    if test_inchi is not None:
        return Chem.MolFromInchi(test_inchi)


def _atom_test(test_inchi=None, test_smiles=None):
    mol = _evaluation_mol_generator(test_inchi=test_inchi, test_smiles=test_smiles)
    return detect_inorganic(mol)


def test_organic():
    """Tests if organic structures are not detected as inorganic.
    Organic atoms:
    hydrogen, carbon, nitrogen, oxygen,
    fluorine, phosphorus, sulfur, chlorine, bromine, iodine
    """
    assert _atom_test(test_inchi="InChI=1S/H") == False
    assert _atom_test(test_inchi="InChI=1S/C") == False
    assert _atom_test(test_inchi="InChI=1S/N") == False
    assert _atom_test(test_inchi="InChI=1S/O") == False
    assert _atom_test(test_inchi="InChI=1S/F") == False
    assert _atom_test(test_inchi="InChI=1S/P") == False
    assert _atom_test(test_inchi="InChI=1S/S") == False
    assert _atom_test(test_inchi="InChI=1S/Cl") == False
    assert _atom_test(test_inchi="InChI=1S/Br") == False
    assert _atom_test(test_inchi="InChI=1S/I") == False


def test_inorganic():
    """Tests if inorganic structures are detected as inorganic.
    Example atoms:
    aluminum, selenium, sodium, magnesium
    """
    assert _atom_test(test_inchi="InChI=1S/Al") == True
    assert _atom_test(test_inchi="InChI=1S/Se") == True
    assert _atom_test(test_inchi="InChI=1S/Na") == True
    assert _atom_test(test_inchi="InChI=1S/Mg") == True


def test_organic_inorganic_combination():
    """Tests if inorganic atoms are detected in combination with organic
    atoms.
    """
    assert _atom_test(test_smiles="c1ccccc1C(=O)O[Ca]OC(=O)c1ccccc1") == True


def test_organic_organic_combination():
    """Tests if a combination of organic atoms doesn't resolve in a
    detection of a inorganic structure.
    Tested here: nitro group
    """
    assert _atom_test(test_smiles="[N](=O)(=O)O") == False
