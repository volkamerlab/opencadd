"""
test for the module `remove_salts`
"""
import pytest
import sys
import rdkit

from rdkit import Chem

from opencadd.compounds.standardization.remove_salts import remove_salts


def _evaluation_mol_generator(test_smiles=None, test_inchi=None):
    """Creates mol files directly with rdkits functions for evaluation."""
    if test_smiles is not None:
        return Chem.MolFromSmiles(test_smiles)
    if test_inchi is not None:
        return Chem.MolFromInchi(test_inchi)


def _molecule_test(test_inchi=None, test_smiles=None):
    return Chem.MolToInchi(
        remove_salts(
            _evaluation_mol_generator(test_inchi=test_inchi, test_smiles=test_smiles)
        )
    )


def test_structure():
    """Only C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-] should be
    left after stripping salts.
    """
    assert (
        _molecule_test(
            test_smiles="C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-].CCCCCCCCCCCCCCCCCC(=O)O.OCC(O)C1OC(=O)C(=C1O)O)"
        )
        == "InChI=1S/C6H10N6O4/c7-2-6(5(13)14,3-12(15)16)1-4-8-10-11-9-4/h1-3,7H2,(H2,8,9,10,11,13,14)/p-1"
    )


def test_single_salts():
    """All salt fragments should be detected and stripped."""
    assert (
        _molecule_test(
            test_smiles="[Al].N.[Ba].[Bi].Br.[Ca].Cl.F.I.[K].[Li].[Mg].[Na].[Ag].[Sr].S.O.[Zn]"
        )
        == ""
    )


def test_complex_salts():
    """Complex salts, contained in salts.tsv should be detected."""
    assert (
        _molecule_test(test_smiles="OC(C(O)C(=O)O)C(=O)O.O=C1NS(=O)(=O)c2ccccc12") == ""
    )


def test_custom_dictionary():
    """Configuration of a custom dictionary, by defining one, should
    work.
    """
    assert (
        Chem.MolToInchi(
            remove_salts(
                _evaluation_mol_generator(test_smiles="[Al].N.[Ba].[Bi]"),
                dictionary=False,
                defnData="[Al]",
            )
        )
        == "InChI=1S/Ba.Bi.H3N.2H/h;;1H3;;"
    )
