"""
This module provides all validation methods included in rdMolStandardize
"""
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.MolStandardize import rdMolStandardize

__all__ = [
    "check_valency",
    "validate_allowed_atoms",
    "validate_disallowed_atoms",
    "validate_default",
    "validate_if_no_atom",
    "validate_fragment",
    "validate_neutrality",
    "validate_isotopes",
]


def check_valency(mol):
    """Validates the valency of every atom in the molecule.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    message: str
    """
    message = rdMolStandardize.RDKitValidation().validate(mol)
    return message or None


def validate_allowed_atoms(mol, atomlist):
    """Validates by a list of allowed atoms.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.
    atomlist: list of int
        The atomic number of the allowed atoms.

    Returns
    -------
    message: str
        Error Message with every atom that is not on atomlist.
    """
    if len(atomlist) == 0:
        return
    else:
        allowed_atoms = [Atom(i) for i in atomlist]
        vm = rdMolStandardize.AllowedAtomsValidation(allowed_atoms)
        message = vm.validate(mol)
        return message or None


def validate_disallowed_atoms(mol, atomlist):
    """
    Validates by a list of disallowed atoms.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.
    atomlist: list of int
        The atomic number of the disallowed atoms.

    Returns
    -------
    message: str
        Error Message with every atom that is on the atomlist.

    """
    if len(atomlist) == 0:
        return
    else:
        disallowed_atoms = [Atom(i) for i in atomlist]
        vm = rdMolStandardize.DisallowedAtomsValidation(disallowed_atoms)
        message = vm.validate(mol)
        return message or None


def validate_default(mol):
    """Performs the default validation which occurs in MolVS.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    message: str

    Notes
    -----
    Default is to do all following validations:
    validate_if_no_atom, validate_fragment, validate_neutrality,
    validate_isotopes
    """
    message = rdMolStandardize.MolVSValidation().validate(mol)
    return message or None


def validate_if_no_atom(mol):
    """Logs an error if the molecule has zero atoms.

    If the molecule has no atoms, no subsequent validations will run.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    message: str
    """
    validation = [rdMolStandardize.NoAtomValidation()]
    message = rdMolStandardize.MolVSValidation(validation).validate(mol)
    return message or None


def validate_fragment(mol):
    """Logs if certain fragments are present.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    message: str
    """
    validation = [rdMolStandardize.FragmentValidation()]
    message = rdMolStandardize.MolVSValidation(validation).validate(mol)
    return message or None


def validate_neutrality(mol):
    """Logs if not an overall neutral system.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    message: str
    """
    validation = [rdMolStandardize.NeutralValidation()]
    message = rdMolStandardize.MolVSValidation(validation).validate(mol)
    return message or None


def validate_isotopes(mol):
    """Logs if molecule contains isotopes.

    Parameters
    ----------
    mol: rdkit.Chem.Mol
        A molecule.

    Returns
    -------
    message: str
    """
    validation = [rdMolStandardize.IsotopeValidation()]
    message = rdMolStandardize.MolVSValidation(validation).validate(mol)
    return message or None
