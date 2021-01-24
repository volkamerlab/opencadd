"""
test for the module `normalize_molecules`

derived from MolVS's tests: https://github.com/mcs07/MolVS/blob/master/tests/test_normalize.py
"""
import pytest

from rdkit import Chem

from opencadd.compounds.standardization import normalize_molecules


def normalization_for_smiles(smiles):
    """Does normalization after converting a SMILES string into a mol."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = normalize_molecules.normalize(mol)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_nitro():
    """Test nitro group normalozation."""
    assert normalization_for_smiles("CN(=O)=O") == "C[N+](=O)[O-]"


def test_sulfoxide():
    """Test sulfoxide normalization."""
    assert normalization_for_smiles("CS(C)=O") == "C[S+](C)[O-]"


def test_sulfone():
    """Test sulfone normalization."""
    assert normalization_for_smiles("C[S+2]([O-])([O-])O") == "CS(=O)(=O)O"


def test_1_3_charge_recombination():
    """Test 1,3-seperated charges are recombined"""
    assert normalization_for_smiles("CC([O-])=[N+](C)C") == "CC(=O)N(C)C"


def test_1_3_charge_recombination_aromatic():
    """Test 1,3-separated charges are recombined."""
    assert normalization_for_smiles("C[n+]1ccccc1[O-]") == "Cn1ccccc1=O"


def test_1_3_charge_recombination_exception():
    """Test a case where 1,3-separated charges should not be recombined."""
    assert (
        normalization_for_smiles("CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]")
        == "CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]"
    )


def test_1_5_charge_recombination():
    """Test 1,5-separated charges are recombined."""
    assert normalization_for_smiles("C[N+](C)=C\\C=C\\[O-]") == "CN(C)C=CC=O"


def test_1_5_charge_recombination_exception():
    """Test a case where 1,5-separated charges should not be recombined."""
    assert (
        normalization_for_smiles("C[N+]1=C2C=[N+]([O-])C=CN2CCC1")
        == "C[N+]1=C2C=[N+]([O-])C=CN2CCC1"
    )
