"""
Tests for opencadd.io.biopython
"""

from pathlib import Path

import pytest
from Bio.PDB.Structure import Structure

from opencadd.io.biopython import PdbToBiopython

PATH_TEST_DATA = Path(__name__).parent / "opencadd" / "tests" / "data" / "io"


def _biopython_format_tests(biopython_structure):
    """
    Default tests: Check if output is an biopython structure object.

    Parameter
    ---------
    biopython_structure : Bio.PDB.Structure.Structure
        Biopython structure object.
    """

    assert isinstance(biopython_structure, Structure)


@pytest.mark.parametrize("pdb_file, n_residues", [(PATH_TEST_DATA / "2itz.pdb", 370)])
def test_from_pdb_file(pdb_file, n_residues):
    """
    Test loading pdb files as biopython structure object.

    Parameters
    ----------
    pdb_file : pathlib.Path or str
        Path to pdb file.
    """

    pdb_to_bpy = PdbToBiopython()
    biopython_structure = pdb_to_bpy.from_file(pdb_file)

    _biopython_format_tests(biopython_structure)
    assert len(list(biopython_structure.get_residues())) == n_residues
