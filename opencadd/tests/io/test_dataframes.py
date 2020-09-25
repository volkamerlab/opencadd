"""
Tests for opencadd.io.dataframes
"""

from pathlib import Path

import pandas as pd
import pytest

from opencadd.io import DataFrame
from opencadd.io.schema import DATAFRAME_COLUMNS

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
            (PATH_TEST_DATA / "2itz.pdb", "mol2"),
            (PATH_TEST_DATA / "2itz.pdb", "unknown format"),
            (PATH_TEST_DATA / "2itz_chainA_protein.mol2", "pdb"),
            (PATH_TEST_DATA / "2itz_chainA_protein.mol2", "unknown format"),
            (PATH_TEST_DATA / "invalid_data", "mol2"),
            (PATH_TEST_DATA / "invalid_data", "pdb"),
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

    @pytest.mark.parametrize("pdb_file", [PATH_TEST_DATA / "2itz.pdb"])
    def test_from_pdb_file(self, pdb_file):
        """
        Test loading pdb files as DataFrame.

        Parameters
        ----------
        pdb_file : pathlib.Path or str
            Path to pdb file.
        """

        df = DataFrame._from_pdb_file(pdb_file)
        self._dataframe_format_tests(df)

    @pytest.mark.parametrize("pdb_file", [PATH_TEST_DATA / "2itz.pdb"])
    def test_from_pdb_text(self, pdb_file):
        """
        Test loading pdb file contents (text) as DataFrame.

        Parameters
        ----------
        pdb_file : pathlib.Path or str
            Path to pdb file.
        """

        # Let's load a file's content as string (text) to simulate example input data
        with open(pdb_file, "r") as f:
            pdb_text = f.read()

        df = DataFrame._from_pdb_text(pdb_text)
        self._dataframe_format_tests(df)

    @pytest.mark.parametrize(
        "mol2_file",
        [PATH_TEST_DATA / "2itz_chainA_protein.mol2", PATH_TEST_DATA / "2itz_protein.mol2"],
    )
    def test_from_mol2_file(self, mol2_file):
        """
        Test loading mol2 files as DataFrame.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
            Path to mol2 file.
        """

        df = DataFrame._from_mol2_file(mol2_file)
        self._dataframe_format_tests(df)

    @pytest.mark.parametrize(
        "mol2_file",
        [PATH_TEST_DATA / "2itz_chainA_protein.mol2", PATH_TEST_DATA / "2itz_protein.mol2"],
    )
    def test_from_mol2_text(self, mol2_file):
        """
        Test loading mol2 file contents (text) as DataFrame.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
            Path to mol2 file.
        """

        # Let's load a file's content as string (text) to simulate example input data
        with open(mol2_file, "r") as f:
            mol2_text = f.read()

        df = DataFrame._from_mol2_text(mol2_text)
        self._dataframe_format_tests(df)

    def _dataframe_format_tests(self, dataframe):
        """
        Default tests: Check if output is a DataFrame, if column names (and their order) and
        the column dtypes are correct.

        Parameter
        ---------
        dataframe : pandas.DataFrame
            DataFrame with structural data.
        """

        assert isinstance(dataframe, pd.DataFrame)

        # Create DataFrame for default dataframe columns
        default_dataframe_columns = pd.DataFrame(
            DATAFRAME_COLUMNS["default"], columns=["name", "dtype"]
        )
        assert dataframe.columns.to_list() == default_dataframe_columns["name"].to_list()
        assert dataframe.dtypes.to_list() == default_dataframe_columns["dtype"].to_list()
