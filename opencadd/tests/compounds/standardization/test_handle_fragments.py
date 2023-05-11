"""
test for the module `handle fragments`

derived from MolVS's tests: https://github.com/mcs07/MolVS/blob/master/tests/test_fragment.py
"""
import pytest
import sys

from rdkit import Chem

from opencadd.compounds.standardization import handle_fragments


def _remove_fragment_smiles(smiles):
    """Utility function that returns the result SMILES after
    remove_fragments is applied to given a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    mol = handle_fragments.remove_fragments(mol)
    return Chem.MolToSmiles(mol)


def test_remove_single_salt():
    """Single Salt removal."""
    assert _remove_fragment_smiles("CN(C)C.Cl") == "CN(C)C"


def test_remove_multiple_salts():
    """Multiple salt removal."""
    assert _remove_fragment_smiles("CN(C)C.Cl.Cl.Br") == "CN(C)C"


def test_fragment_patterns():
    """FragmentPatterns should match entire fragments only, matches
    within larger fragments should be left.
    """
    assert _remove_fragment_smiles("CN(Br)Cl") == "CN(Cl)Br"
    assert _remove_fragment_smiles("CN(Br)Cl.Cl") == "CN(Cl)Br"


def test_charged_salts():
    """Charged salts."""
    assert _remove_fragment_smiles("C[NH+](C)(C).[Cl-]") == "C[NH+](C)C"


def test_last_match():
    """Last match should be left."""
    assert _remove_fragment_smiles("CC(=O)O.[Na]") == "CC(=O)O"


def test_left_identical():
    """Multiple identical last fragments should all be left."""
    assert _remove_fragment_smiles("Br.Br") == "Br.Br"


def test_remove_multiple_fragment():
    """Test multiple fragment removal."""
    assert (
        _remove_fragment_smiles("[Na+].OC(=O)Cc1ccc(CN)cc1.OS(=O)(=O)C(F)(F)F")
        == "NCc1ccc(CC(=O)O)cc1"
    )


def test_1_4_Dioxiane():
    """1,4-Dioxane should be removed."""
    assert _remove_fragment_smiles("c1ccccc1O.O1CCOCC1") == "Oc1ccccc1"


def test_benzene():
    """Benzene should be removed."""
    assert _remove_fragment_smiles("c1ccccc1.CCCBr") == "CCCBr"


def test_remove_various_fragments():
    """Various fragments should be removed."""
    assert (
        _remove_fragment_smiles(
            "CC(NC1=CC=C(O)C=C1)=O.CCCCC.O.CCO.CCCO.C1CCCCC1.C1CCCCCC1"
        )
        == "CC(=O)Nc1ccc(O)cc1"
    )
