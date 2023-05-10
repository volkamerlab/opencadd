"""
Tests for opencadd.io.biopython
"""

from pathlib import Path

import pytest
from Bio.PDB.Structure import Structure

from opencadd.io.biopython import Biopython

PATH_TEST_DATA = Path(__name__).parent / "opencadd" / "tests" / "data" / "io"


class TestsBiopython:
    """
    Tests the DataFrame class methods.
    """

    @pytest.mark.parametrize(
        "filepath",
        [PATH_TEST_DATA / "2itz.pdb"],
    )
    def test_from_file(self, filepath):
        """
        Test if input produces a biopython structure object.
        """

        biopython_structure = Biopython.from_file(filepath)
        self._biopython_format_tests(biopython_structure)

    @pytest.mark.parametrize(
        "filepath",
        [PATH_TEST_DATA / "2itz_chainA_protein.mol2"],
    )
    def test_from_file_raises(self, filepath):
        """
        Test if input produces a ValueError for invalid inputs.
        """

        with pytest.raises(ValueError):
            Biopython.from_file(filepath)

    @pytest.mark.parametrize("pdb_file, n_residues", [(PATH_TEST_DATA / "2itz.pdb", 370)])
    def test_from_pdb_file(self, pdb_file, n_residues):
        """
        Test loading pdb files as biopython structure object.

        Parameters
        ----------
        pdb_file : pathlib.Path or str
            Path to pdb file.
        """

        biopython_structure = Biopython._from_pdb_file(pdb_file)

        self._biopython_format_tests(biopython_structure)
        assert len(list(biopython_structure.get_residues())) == n_residues

    def _biopython_format_tests(self, biopython_structure):
        """
        Default tests: Check if output is an biopython structure object.

        Parameter
        ---------
        biopython_structure : Bio.PDB.Structure.Structure
            Biopython structure object.
        """

        assert isinstance(biopython_structure, Structure)
