"""
Tests for opencadd.io.rdkit
"""

from pathlib import Path

import pytest
from rdkit import Chem

from opencadd.io import Rdkit

PATH_TEST_DATA = Path(__name__).parent / "opencadd" / "tests" / "data" / "io"


class TestsRdkitMol:
    """
    Tests the RdkitMol class methods.
    """

    @pytest.mark.parametrize(
        "filepath, compute2d",
        [
            (PATH_TEST_DATA / "2itz_chainA_ligand.mol2", True),
            (PATH_TEST_DATA / "2itz_chainA_ligand.mol2", False),
        ],
    )
    def test_from_file(self, filepath, compute2d):
        """
        Test if input produces an RDKit molecule.
        """

        rmol = Rdkit.from_file(filepath, compute2d)
        self._rdkit_format_tests(rmol)

    @pytest.mark.parametrize(
        "filepath, format, compute2d",
        [
            (PATH_TEST_DATA / "2itz_chainA_ligand.mol2", "mol2", True),
            (PATH_TEST_DATA / "2itz_chainA_ligand.mol2", "mol2", False),
        ],
    )
    def test_from_text(self, filepath, format, compute2d):
        """
        Test if input produces an RDKit molecule.
        Note: Use file as test function input; file content will be read as string,
        which is the input for the class method to be tested here!
        """

        # Let's load a file's content as string (text) to simulate example input data
        with open(filepath, "r") as f:
            text = f.read()

        rmol = Rdkit.from_text(text, format, compute2d)
        self._rdkit_format_tests(rmol)

    @pytest.mark.parametrize(
        "mol2_file",
        [PATH_TEST_DATA / "2itz_chainA_ligand.mol2"],
    )
    def test_from_mol2_file(self, mol2_file):
        """
        Test loading mol2 files as RDKit molecule.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
            Path to mol2 file.
        """

        rmol = Rdkit._from_mol2_file(mol2_file)
        self._rdkit_format_tests(rmol)

    @pytest.mark.parametrize(
        "mol2_file",
        [PATH_TEST_DATA / "2itz_chainA_ligand.mol2"],
    )
    def test_from_mol2_text(self, mol2_file):
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

        rmol = Rdkit._from_mol2_text(mol2_text)
        self._rdkit_format_tests(rmol)

    def _rdkit_format_tests(self, rmol):
        """
        Default tests: Check if output is an RDKit molecule.

        Parameter
        ---------
        rmol : Chem.rdchem.Mol
            RDKit molecule.
        """

        assert isinstance(rmol, Chem.rdchem.Mol)
