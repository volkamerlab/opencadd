"""
test for the module `handle_charges`

derived from MolVS's tests: https://github.com/mcs07/MolVS/blob/master/tests/test_charge.py
"""
import pytest
import sys

from rdkit import Chem

from opencadd.compounds.standardization import handle_charges


def _uncharge_smiles(smiles):
    """Utility function that returns the uncharged SMILES for a given
    SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = handle_charges.uncharge(mol)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_neutralization():
    """Test neutralization of ionized acids and bases."""
    assert (
        _uncharge_smiles("C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])")
        == "NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O"
    )


def test_zwitterion():
    """Test preservation of zwitterion."""
    assert _uncharge_smiles("n(C)1cc[n+]2cccc([O-])c12") == "Cn1cc[n+]2cccc([O-])c12"


def test_choline():
    """Choline should be left with a positive charge."""
    assert _uncharge_smiles("C[N+](C)(C)CCO") == "C[N+](C)(C)CCO"


def test_hydrogen():
    """This should have the hydrogen removed to give deanol as a charge parent."""
    assert _uncharge_smiles("C[NH+](C)CCO") == "CN(C)CCO"


def test_neutrality():
    """Overall system is already neutral."""
    assert _uncharge_smiles("[Na+].O=C([O-])c1ccccc1") == "O=C([O-])c1ccccc1.[Na+]"


def test_benzoate():
    """Benzoate ion to benzoic acid."""
    assert _uncharge_smiles("O=C([O-])c1ccccc1") == "O=C(O)c1ccccc1"


def test_histidine():
    """Charges in histidine should be neutralized."""
    assert _uncharge_smiles("[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]") == "NC(Cc1cnc[nH]1)C(=O)O"


def test_fragment_neutralization():
    """Neutralize both fragments."""
    assert _uncharge_smiles("C[NH+](C)(C).[Cl-]") == "CN(C)C.Cl"


def test_oxigen_neutralisation():
    """Neutralise one oxygen."""
    assert _uncharge_smiles("[N+](=O)([O-])[O-]") == "O=[N+]([O-])[O-]"


def test_prefer_organic_fragments():
    """Smaller organic fragment should be chosen over larger inorganic fragment."""
    assert _uncharge_smiles("[N+](=O)([O-])[O-].[CH2]") == "O=[N+]([O-])[O-].[CH2]"


def test_oxygen_balancing():
    """Single oxygen should be protonated, the other left to balance the positive nitrogen."""
    assert _uncharge_smiles("C[N+](C)(C)CC([O-])C[O-]") == "C[N+](C)(C)CC([O-])CO"


def test_strongest_acid():
    """Strongest acid should be left ionized."""
    assert (
        _uncharge_smiles("[O-]C(=O)C[n+]1ccn2cccc([O-])c12")
        == "O=C([O-])C[n+]1ccn2cccc(O)c21"
    )


def test_charge_neutralization():
    """All charges should be neutralized."""
    assert _uncharge_smiles("[NH+](C)(C)CC([O-])C[O-]") == "CN(C)CC(O)CO"


def test_uncharge():
    """All charges should be neutralized."""
    assert _uncharge_smiles("CNCC([O-])C[O-]") == "CNCC(O)CO"


# Tests for Reionize


def _reionize_smiles(smiles):
    """Utility function that returns the uncharged SMILES for a given
    SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    mol = handle_charges.reionize(mol)
    if mol:
        return Chem.MolToSmiles(mol)


def test_proton_to_weak_acid():
    """Test reionizer moves proton to weaker acid."""
    assert (
        _reionize_smiles("C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O")
        == "O=S(O)c1ccc(S(=O)(=O)[O-])cc1"
    )


def test_charged_carbon():
    """Test charged carbon doesn't get recognised as
    alpha-carbon-hydrogen-keto.
    """
    assert _reionize_smiles("CCOC(=O)C(=O)[CH-]C#N") == "CCOC(=O)C(=O)[CH-]C#N"
