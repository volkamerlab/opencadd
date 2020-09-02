"""
opencadd.io.dataframes

Defines classes that convert different input types into DataFrames.
"""

from pathlib import Path

from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd

from .schema import DATAFRAME_COLUMNS, PDB_COLUMNS, MOL2_COLUMNS


class ToDataFrame:
    """
    Base class for DataFrame generation from different input formats.
    """

    def _format_dataframe(self, dataframe, verbose=False):
        """
        Streamline the output DataFrame's column name order and dtypes accross all input formats.

        Parameter
        ---------
        dataframe : pandas.DataFrame
            DataFrame with structural data.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Formatted DataFrame with structural data.
        """

        if verbose:
            dataframe_columns = DATAFRAME_COLUMNS["default"] + DATAFRAME_COLUMNS["verbose"]
        else:
            dataframe_columns = DATAFRAME_COLUMNS["default"]

        # If default column names are missing, add empty columns
        for (column_name, column_dtype) in dataframe_columns:

            if column_name not in dataframe.columns:

                if column_dtype == "object":
                    dataframe.insert(len(dataframe.columns), column_name, None)
                elif column_dtype == "float64":
                    dataframe.insert(len(dataframe.columns), column_name, np.nan)
                elif column_dtype == "Int64":
                    dataframe.insert(len(dataframe.columns), column_name, pd.NA)
                else:
                    raise KeyError(f"Column dtype {column_dtype} is not implemented.")

        # Set default dtypes
        column_dtypes_dict = {
            column_name: column_dtype for (column_name, column_dtype) in dataframe_columns
        }
        dataframe = dataframe.astype(column_dtypes_dict)

        # Set default columns order
        column_names_list = [column_name for (column_name, column_dtype) in dataframe_columns]
        dataframe = dataframe[column_names_list]

        return dataframe.dropna(axis=1)


class PdbToDataFrame(ToDataFrame):
    """
    Parse a pdb file or pdb text into a DataFrame.
    """

    def from_file(self, pdb_file, verbose=False):
        """
        Get structural data from pdb file.

        Parameters
        ----------
        pdb_file : pathlib.Path or str
            Path to pdb file.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        pdb_file = Path(pdb_file)

        if not pdb_file.exists():
            raise FileNotFoundError(f"File {pdb_file} does not exist.")

        with open(pdb_file, "r") as f:
            text = f.read()

        return self.from_text(text, verbose)

    def from_text(self, pdb_text, verbose=False):
        """
        Get structural data from pdb text.

        Parameters
        ----------
        pdb_text : str
            Pdb file content from KLIFS database.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        dict of pandas.DataFrame
            Structural data
        """

        # Use biopandas to parse the pdb format and return DataFrames
        ppdb = PandasPdb()
        pdb_dict = ppdb._construct_df(pdb_text.splitlines(True))

        # Concatenate ATOM and HETATM entries
        pdb_df = pd.concat([pdb_dict["ATOM"], pdb_dict["HETATM"]])

        # Select only columns of interest and rename columns
        pdb_df = pdb_df.iloc[:, list(PDB_COLUMNS.keys())]
        pdb_df.columns = [i[0] for i in PDB_COLUMNS.values()]

        # Format DataFrame
        pdb_df = self._format_dataframe(pdb_df, verbose)

        return pdb_df


class Mol2ToDataFrame(ToDataFrame):
    """
    Parse a mol2 file or mol2 text into a DataFrame.
    """

    def from_file(self, mol2_file, verbose=False):
        """
        Get structural data from mol2 file.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
            Path to mol2 file.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        mol2_file = Path(mol2_file)

        if not mol2_file.exists():
            raise FileNotFoundError(f"File {mol2_file} does not exist.")

        with open(mol2_file, "r") as f:
            text = f.read()

        return self.from_text(text, verbose)

    def from_text(self, mol2_text, verbose=False):
        """
        Get structural data from mol2 text.

        Parameters
        ----------
        mol2_text : str
            Mol2 file content from KLIFS database.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        # Use biopandas to parse the mol2 format and return a DataFrame
        pmol = PandasMol2()
        try:
            mol2_df = pmol._construct_df(
                mol2_text.splitlines(True),
                col_names=[i[0] for i in MOL2_COLUMNS["n_cols_10"].values()],
                col_types=[i[1] for i in MOL2_COLUMNS["n_cols_10"].values()],
            )
        except ValueError:
            mol2_df = pmol._construct_df(
                mol2_text.splitlines(True),
                col_names=[i[0] for i in MOL2_COLUMNS["n_cols_9"].values()],
                col_types=[i[1] for i in MOL2_COLUMNS["n_cols_9"].values()],
            )

        # Infer residue PDB ID and name from substructure name
        mol2_df = self._split_mol2_subst_names(mol2_df)

        # Format DataFrame
        mol2_df = self._format_dataframe(mol2_df, verbose)

        return mol2_df

    def _split_mol2_subst_names(self, mol2_df):
        """
        Split "residue.subst_name" column values from mol2 file and add as 
        "residue.name" and "residue.pdb_id" columns.

        Parameters
        ----------
        mol2_df : pandas.DataFrame
            Structural data.

        Returns
        -------
        pandas.DataFrame
            Structural data, including additional columns for residue name and PDB ID.
        """

        result = mol2_df.apply(
            lambda x: self._split_mol2_subst_name(x["residue.subst_name"], x["atom.type"]), axis=1,
        )
        res_names = [res_name for res_name, res_id in result]
        res_ids = [res_id for res_name, res_id in result]

        mol2_df["residue.name"] = res_names
        mol2_df["residue.pdb_id"] = res_ids

        return mol2_df

    @staticmethod
    def _split_mol2_subst_name(subst_name, atom_type):
        """
        Split "residue.subst_name" values from mol2 file.

        Parameters
        ----------
        subst_name : str
            Substructure name in mol2 file. 
        atom_type : str
            Atom type.

        Returns
        -------
        tuple (str, int)
            Residue name and residue PDB ID.

        Notes
        -----
        Regular example: "ALA1" ("residue.subst_name") will be split to 
        - "ALA" ("residue.name") and 
        - 1 ("residue.pdb_id").
        """

        # Handle "residues" that are elements such as CA or MG.
        if subst_name[:2] == atom_type.upper():
            res_name = subst_name[:2]
            res_id = subst_name[2:]

        # These are amino acid, linkers, compounds, ...
        else:
            res_name = subst_name[:3]
            res_id = subst_name[3:]

        return res_name, res_id
