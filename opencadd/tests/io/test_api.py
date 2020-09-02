"""
Tests for opencadd.io.core
"""

from pathlib import Path

import pytest
import pandas as pd
from rdkit import Chem

from opencadd.io import DataFrame, RdkitMol

PATH_TEST_DATA = Path(__name__).parent / "opencadd" / "tests" / "data" / "io"


class TestsDataFrame:
    """
    Tests the DataFrame class methods.
    """

    @pytest.mark.parametrize(
        "filepath, verbose",
        [
            (PATH_TEST_DATA / "2itz.pdb", False),
            (PATH_TEST_DATA / "2itz.pdb", True),
            (PATH_TEST_DATA / "2itz_chainA_protein.mol2", False),
            (PATH_TEST_DATA / "2itz_protein.mol2", True),
            (PATH_TEST_DATA / "2itz_protein.mol2", False),
            (PATH_TEST_DATA / "2itz_protein.mol2", True),
        ],
    )
    def test_from_file(self, filepath, verbose):
        """
        Test if input produces a DataFrame.
        """

        dataframe = DataFrame.from_file(filepath, verbose)
        isinstance(dataframe, pd.DataFrame)

    @pytest.mark.parametrize(
        "filepath, format, verbose",
        [
            (PATH_TEST_DATA / "2itz.pdb", "pdb", False),
            (PATH_TEST_DATA / "2itz.pdb", "pdb", True),
            (PATH_TEST_DATA / "2itz_chainA_protein.mol2", "mol2", False),
            (PATH_TEST_DATA / "2itz_protein.mol2", "mol2", True),
            (PATH_TEST_DATA / "2itz_protein.mol2", "mol2", False),
            (PATH_TEST_DATA / "2itz_protein.mol2", "mol2", True),
        ],
    )
    def test_from_text(self, filepath, format, verbose):
        """
        Test if input produces a DataFrame. 
        Note: Use file as test function input; file content will be read as string, 
        which is the input for the class method to be tested here!
        """

        # Let's load a file's content as string (text) to simulate example input data
        with open(filepath, "r") as f:
            text = f.read()

        dataframe = DataFrame.from_text(text, format, verbose)
        isinstance(dataframe, pd.DataFrame)

    @pytest.mark.parametrize(
        "filepath, format",
        [
            (PATH_TEST_DATA / "2itz.pdb", "mol2"),  # TODO
            (PATH_TEST_DATA / "2itz.pdb", "unknown format"),
            (PATH_TEST_DATA / "2itz_chainA_protein.mol2", "pdb"),  # TODO
            (PATH_TEST_DATA / "2itz_chainA_protein.mol2", "unknown format"),
            (PATH_TEST_DATA / "invalid_data", "mol2"),  # TODO
            (PATH_TEST_DATA / "invalid_data", "pdb"),  # TODO
            (PATH_TEST_DATA / "invalid_data", "unknown format"),
        ],
    )
    def test_from_text_raises(self, filepath, format):
        """
        Test if input produces a ValueError for invalid inputs. 
        Note: Use file as test function input; file content will be read as string, 
        which is the input for the class method to be tested here!
        """

        # Let's load a file's content as string (text) to simulate example input data
        with open(filepath, "r") as f:
            text = f.read()

        with pytest.raises(ValueError):
            DataFrame.from_text(text, format)


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

        rmol = RdkitMol.from_file(filepath, compute2d)
        isinstance(rmol, Chem.rdchem.Mol)

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

        rmol = RdkitMol.from_text(text, format, compute2d)
        isinstance(rmol, Chem.rdchem.Mol)
