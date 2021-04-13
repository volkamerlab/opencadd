"""
test for the module `handle_hydrogens`
"""
import pytest

from rdkit import Chem

from opencadd.compounds.standardization import handle_hydrogens


def _evaluation_mol_generator(test_smiles=None, test_inchi=None):
    """Creates mol files directly with rdkits functions for evaluation."""
    if test_smiles is not None:
        return Chem.MolFromSmiles(test_smiles, sanitize=False)
    if test_inchi is not None:
        return Chem.MolFromInchi(test_inchi)


def test_remove_hydrogen():
    assert (
        Chem.MolToInchi(
            handle_hydrogens.remove_hydrogens(
                _evaluation_mol_generator(
                    test_inchi="InChI=1S/2C7H6O2.Ca/c2*8-7(9)6-4-2-1-3-5-6;/h2*1-5H,(H,8,9);/q;;+2/p-2"
                )
            )
        )
        == "InChI=1S/2C7H6O2.Ca/c2*8-7(9)6-4-2-1-3-5-6;/h2*1-5H,(H,8,9);/q;;+2/p-2"
    )
