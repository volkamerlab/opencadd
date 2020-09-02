"""
Tests for opencadd.io.dataframes
"""

from pathlib import Path

import pandas as pd
import pytest

from opencadd.io.dataframes import Mol2ToDataFrame, PdbToDataFrame
from opencadd.io.schema import DATAFRAME_COLUMNS

PATH_TEST_DATA = Path(__name__).parent / "opencadd" / "tests" / "data" / "io"


def _dataframe_format_tests(dataframe):
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


@pytest.mark.parametrize("pdb_file", [PATH_TEST_DATA / "2itz.pdb"])
def test_from_pdb_file(pdb_file):
    """
    Test loading pdb files as DataFrame.

    Parameters
    ----------
    pdb_file : pathlib.Path or str
        Path to pdb file.
    """

    pdb_to_df = PdbToDataFrame()
    df = pdb_to_df.from_file(pdb_file)

    _dataframe_format_tests(df)


@pytest.mark.parametrize("pdb_file", [PATH_TEST_DATA / "2itz.pdb"])
def test_from_pdb_text(pdb_file):
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

    pdb_to_df = PdbToDataFrame()
    df = pdb_to_df.from_text(pdb_text)

    _dataframe_format_tests(df)


@pytest.mark.parametrize(
    "mol2_file",
    [PATH_TEST_DATA / "2itz_chainA_protein.mol2", PATH_TEST_DATA / "2itz_protein.mol2"],
)
def test_from_mol2_file(mol2_file):
    """
    Test loading mol2 files as DataFrame.

    Parameters
    ----------
    mol2_file : pathlib.Path or str
        Path to mol2 file.
    """

    mol2_to_df = Mol2ToDataFrame()
    df = mol2_to_df.from_file(mol2_file)

    _dataframe_format_tests(df)


@pytest.mark.parametrize(
    "mol2_file",
    [PATH_TEST_DATA / "2itz_chainA_protein.mol2", PATH_TEST_DATA / "2itz_protein.mol2"],
)
def test_from_mol2_text(mol2_file):
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

    mol2_to_df = Mol2ToDataFrame()
    df = mol2_to_df.from_text(mol2_text)

    _dataframe_format_tests(df)
