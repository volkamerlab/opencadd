"""
Tests for opencadd.io.rdkit
"""

from pathlib import Path

import pytest
from rdkit import Chem

from opencadd.io.rdkit import Mol2ToRdkitMol

PATH_TEST_DATA = Path(__name__).parent / "opencadd" / "tests" / "data" / "io"


def _rdkit_format_tests(rmol):
    """
    Default tests: Check if output is an RDKit molecule.

    Parameter
    ---------
    rmol : Chem.rdchem.Mol
        RDKit molecule.
    """

    assert isinstance(rmol, Chem.rdchem.Mol)


@pytest.mark.parametrize(
    "mol2_file",
    [PATH_TEST_DATA / "2itz_chainA_ligand.mol2"],
)
def test_from_mol2_file(mol2_file):
    """
    Test loading mol2 files as RDKit molecule.

    Parameters
    ----------
    mol2_file : pathlib.Path or str
        Path to mol2 file.
    """

    mol2_to_rmol = Mol2ToRdkitMol()
    rmol = mol2_to_rmol.from_file(mol2_file)

    _rdkit_format_tests(rmol)


@pytest.mark.parametrize(
    "mol2_file",
    [PATH_TEST_DATA / "2itz_chainA_ligand.mol2"],
)
def test_from_mol2_text(mol2_file):
    """
    Test loading mol2 file contents (text) as RDKit molecule.

    Parameters
    ----------
    mol2_file : pathlib.Path or str
        Path to mol2 file.
    """

    # Let's load a file's content as string (text) to simulate example input data
    with open(mol2_file, "r") as f:
        mol2_text = f.read()

    mol2_to_rmol = Mol2ToRdkitMol()
    rmol = mol2_to_rmol.from_text(mol2_text)

    _rdkit_format_tests(rmol)
