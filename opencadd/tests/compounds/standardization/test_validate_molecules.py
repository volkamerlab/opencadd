"""
test for the module `validate_molecules`

Partially derived from MolVS's tests and RDKIT MolStandardize tutorial:
https://github.com/mcs07/MolVS/blob/master/tests/test_validate.py
https://github.com/susanhleung/rdkit/blob/dev/GSOC2018_MolVS_Integration/rdkit/Chem/MolStandardize/tutorial/MolStandardize.ipynb
"""
import pytest

from rdkit import Chem

from opencadd.compounds.standardization import validate_molecules


def _evaluation_mol_generator(test_smiles=None, test_inchi=None):
    """Creates mol files directly with rdkits functions for evaluation."""
    if test_smiles is not None:
        return Chem.MolFromSmiles(test_smiles, sanitize=False)
    if test_inchi is not None:
        return Chem.MolFromInchi(test_inchi)


def test_no_atom():
    """NoAtomValidation should log due to the lack of any atoms."""
    assert validate_molecules.validate_default(
        _evaluation_mol_generator(test_smiles="")
    ) == ["ERROR: [NoAtomValidation] Molecule has no atoms"]


def test_fragment_dichloroethane():
    """FragmentValidation should identify 1,2-dichloroethane."""
    assert validate_molecules.validate_fragment(
        _evaluation_mol_generator(test_smiles="ClCCCl.c1ccccc1O")
    ) == ["INFO: [FragmentValidation] 1,2-dichloroethane is present"]


def test_fragment_dimethoxyethane():
    """FragmentValidation should identify 1,2-dimethoxyethane."""
    assert validate_molecules.validate_fragment(
        _evaluation_mol_generator(test_smiles="COCCOC.CCCBr")
    ) == ["INFO: [FragmentValidation] 1,2-dimethoxyethane is present"]


def test_neutrality():
    """NeutralValidation should identify net overall charge."""
    assert validate_molecules.validate_neutrality(
        _evaluation_mol_generator(test_smiles="O=C([O-])c1ccccc1")
    ) == ["INFO: [NeutralValidation] Not an overall neutral system (-1)"]
    assert validate_molecules.validate_neutrality(
        _evaluation_mol_generator(test_smiles="CN=[NH+]CN=N")
    ) == ["INFO: [NeutralValidation] Not an overall neutral system (+1)"]


def test_isotope():
    """IsotopeValidation should identify atoms with isotope labels."""
    assert validate_molecules.validate_isotopes(
        _evaluation_mol_generator(test_smiles="[13CH4]")
    ) == ["INFO: [IsotopeValidation] Molecule contains isotope 13C"]
    assert validate_molecules.validate_isotopes(
        _evaluation_mol_generator(test_smiles="[2H]C(Cl)(Cl)Cl")
    ) == ["INFO: [IsotopeValidation] Molecule contains isotope 2H"]
    assert validate_molecules.validate_isotopes(
        _evaluation_mol_generator(test_smiles="[2H]OC([2H])([2H])[2H]")
    ) == ["INFO: [IsotopeValidation] Molecule contains isotope 2H"]


def test_valency():
    """check_valency should validate the valency of every atom in the
    input molecule.
    """
    assert validate_molecules.check_valency(
        _evaluation_mol_generator(test_smiles="CO(C)C")
    ) == [
        "INFO: [ValenceValidation] Explicit valence for atom # 1 O, 3, is greater than permitted"
    ]


def test_validate_allowed_atoms():
    """validate_allowed_atoms should accept as input a list of atoms,
    anything not on the list should throw an error.
    """
    assert validate_molecules.validate_allowed_atoms(
        _evaluation_mol_generator(test_smiles="CC(=O)CF"), atomlist=[6, 7, 8]
    ) == ["INFO: [AllowedAtomsValidation] Atom F is not in allowedAtoms list"]


def test_validate_disallowed_atoms():
    """validate_allowed_atoms should accept as input a list of atoms,
    anything not on the list should throw an error.
    """
    assert validate_molecules.validate_disallowed_atoms(
        _evaluation_mol_generator(test_smiles="CC(=O)CF"), atomlist=[9, 17, 35]
    ) == ["INFO: [DisallowedAtomsValidation] Atom F is in disallowedAtoms list"]
