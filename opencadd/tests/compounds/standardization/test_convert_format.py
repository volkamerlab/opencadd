"""
test for the module `convert_format`
"""
import pytest
import sys
import rdkit
import os
from pathlib import Path
from rdkit import Chem

from opencadd.compounds.standardization import convert_format


def _evaluation_mol_generator(test_smiles=None, test_inchi=None):
    """Creates mol files directly with rdkits functions for evaluation."""
    if test_smiles is not None:
        return Chem.MolFromSmiles(test_smiles)
    if test_inchi is not None:
        return Chem.MolFromInchi(test_inchi)


def _evaluation_inchi(test_result):
    test_result = Chem.MolFromSmiles(test_result)
    test_result = Chem.MolToInchi(test_result)
    return test_result


def _test_path(fn):
    """Leads to files saved in the data folder

    Parameters
    ----------
    fn: str
        The whole filename.

    Returns
    -------
    The path of the file in the current working system.
    """
    return Path(__file__).parent / "data" / fn


def test_convert_smiles_to_mol(test_smiles="C(C1C(C(C(C(O1)O)O)O)O)O"):
    """Tests if the created file is a mol file."""
    test_result = convert_format.convert_smiles_to_mol(test_smiles)
    assert isinstance(test_result, rdkit.Chem.rdchem.Mol) == True


def test_convert_inchi_to_mol(
    test_inchi="InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1",
):
    """Tests if the created file is a mol file."""
    test_result = convert_format.convert_inchi_to_mol(test_inchi)
    assert isinstance(test_result, rdkit.Chem.rdchem.Mol) == True


def test_convert_mol_to_smiles(
    test_inchi="InChI=1S/C8H7O4S.Na/c9-8(10)7-3-1-6(2-4-7)5-13(11)12;/h1-4H,5H2,(H,9,10);/q;+1/p-1",
):
    """Tests if the created file matches the file it originated from."""
    test_result = convert_format.convert_mol_to_smiles(
        _evaluation_mol_generator(test_inchi=test_inchi), canonical=False
    )
    test_result = _evaluation_inchi(test_result)
    assert test_result == test_inchi


def test_convert_mol_to_inchi(
    test_inchi="InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1",
):
    """Tests if the created file matches the file it originated from."""
    test_result = convert_format.convert_mol_to_inchi(
        _evaluation_mol_generator(test_inchi=test_inchi)
    )
    assert test_result == test_inchi


def test_convert_sdf_to_mol_array():
    fn = "new_mol"
    mol_array = convert_format.convert_sdf_to_mol_array(str(_test_path(fn)))
    mol = mol_array[0]
    assert isinstance(mol, rdkit.Chem.rdchem.Mol) == True


def test_convert_mol_to_sdf():
    fn = "new_mol"
    mol_array = convert_format.convert_sdf_to_mol_array(str(_test_path(fn)))
    convert_format.convert_mol_to_sdf(
        mol_array, fn=str(_test_path("result_of_test_convert_mol_to_sdf.sdf"))
    )
    assert (
        os.path.isfile(str(_test_path("result_of_test_convert_mol_to_sdf.sdf"))) == True
    )
