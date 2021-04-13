"""
test for the module `disconnect_metals`

derived from MolVS's tests:
https://github.com/mcs07/MolVS/blob/master/tests/test_metal.py
"""
import pytest
import sys

from rdkit import Chem

from opencadd.compounds.standardization import disconnect_metals


def _disconnect_metals_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = disconnect_metals(mol)
    if mol:
        return Chem.MolToSmiles(mol)


def test_disconnect_metals1():
    assert (
        _disconnect_metals_smiles("NC(CC(=O)O)C(=O)[O-].O.O.[Na+]")
        == "NC(CC(=O)O)C(=O)[O-].O.O.[Na+]"
    )


def test_covalent_metal():
    """Test if covalent metal is disconnected."""
    assert _disconnect_metals_smiles("CCC(=O)O[Na]") == "CCC(=O)[O-].[Na+]"


def test_no_accidental_deletion():
    """Test metal ion is untouched."""
    assert _disconnect_metals_smiles("CCC(=O)[O-].[Na+]") == "CCC(=O)[O-].[Na+]"


def test_dimethylmercury():
    """Test dimethylmercury is not disconnected."""
    assert _disconnect_metals_smiles("C[Hg]C") == "C[Hg]C"


def test_zirconium():
    """Test zirconium (IV) ethoxide."""
    assert (
        _disconnect_metals_smiles("CCO[Zr](OCC)(OCC)OCC")
        == "CC[O-].CC[O-].CC[O-].CC[O-].[Zr+4]"
    )
