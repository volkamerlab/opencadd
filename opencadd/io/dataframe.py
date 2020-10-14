"""
opencadd.io.dataframes

Defines classes that convert structural data into DataFrames.
"""

import logging

from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import numpy as np
import pandas as pd

from .core import _Base
from .schema import DATAFRAME_COLUMNS, PDB_COLUMNS, MOL2_COLUMNS

logger = logging.getLogger(__name__)


class DataFrame(_Base):
    """
    Parse a structure as a DataFrame.
    """

    @classmethod
    def from_file(cls, filepath, verbose=False, **kwargs):
        """
        Load structures as DataFrame from file.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to structure file: pdb and mol2 files only.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structure as DataFrame.
        """

        filepath = cls._convert_filepath(filepath)

        if filepath.suffix == ".mol2":
            return cls._from_mol2_file(filepath, verbose)
        elif filepath.suffix == ".pdb":
            return cls._from_pdb_file(filepath, verbose)
        else:
            raise ValueError(f"The {filepath.suffix} format is not supported or invalid.")

        return dataframe

    @classmethod
    def from_text(cls, text, ext, verbose=False, **kwargs):
        """
        Load structures as DataFrame from text (file content as string).

        Parameters
        ----------
        text : str
            Structure file content as string.
        ext : str
            Structure format: "pdb" or "mol2".
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structure as DataFrame.
        """

        if ext == "mol2":
            return cls._from_mol2_text(text, verbose)
        elif ext == "pdb":
            return cls._from_pdb_text(text, verbose)
        else:
            raise ValueError(f"The {ext} format is not supported or invalid.")

    @classmethod
    def _from_pdb_file(cls, pdb_file, verbose=False):
        """
        Get structural data from pdb file as DataFrame.

        Parameters
        ----------
        pdb_file : pathlib.Path
            Path to pdb file.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        text = cls._file_to_text(pdb_file)
        return cls._from_pdb_text(text, verbose)

    @classmethod
    def _from_pdb_text(cls, pdb_text, verbose=False):
        """
        Get structural data from pdb text as DataFrame.

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

        # Set pdb columns (index, name, dtype) as DataFrame
        pdb_columns = pd.DataFrame.from_dict(
            PDB_COLUMNS, orient="index", columns=["name", "dtype"]
        )

        # Use biopandas to parse the pdb format and return DataFrames
        # TODO in the future: BioPandas: wait for pdb equivalent of PandasMol2.read_mol2_from_list
        ppdb = PandasPdb()
        pdb_dict = ppdb._construct_df(pdb_text.splitlines(True))

        # Concatenate ATOM and HETATM entries
        pdb_df = pd.concat([pdb_dict["ATOM"], pdb_dict["HETATM"]]).reset_index(drop=True)

        # Select only columns of interest and rename columns
        pdb_df = pdb_df.iloc[:, pdb_columns.index.to_list()]
        pdb_df.columns = pdb_columns["name"].to_list()

        # Format DataFrame
        pdb_df = cls._format_dataframe(pdb_df, verbose)

        if len(pdb_df) == 0:
            raise ValueError(
                f"No structural data could be loaded. Is the input text in pdb format?"
            )

        return pdb_df

    @classmethod
    def _from_mol2_file(cls, mol2_file, verbose=False):
        """
        Get structural data from mol2 file as DataFrame.

        Parameters
        ----------
        mol2_file : pathlib.Path
            Path to mol2 file.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        text = cls._file_to_text(mol2_file)
        return cls._from_mol2_text(text, verbose)

    @classmethod
    def _from_mol2_text(cls, mol2_text, verbose=False):
        """
        Get structural data from mol2 text as DataFrame.

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

        mol2_text = mol2_text.split("\n")

        # Use biopandas to parse the mol2 format and return a DataFrame
        try:
            pmol = PandasMol2()
            try:
                mol2_df = pmol.read_mol2_from_list(
                    mol2_text, "mol", columns=MOL2_COLUMNS["n_cols_10"]
                ).df
            except ValueError as e:
                if str(e) == "10 columns passed, passed data had 9 columns":
                    mol2_df = pmol.read_mol2_from_list(
                        mol2_text, "mol", columns=MOL2_COLUMNS["n_cols_9"]
                    ).df
                else:
                    raise e
        except UnboundLocalError as e:
            if str(e) == "local variable 'first_idx' referenced before assignment":
                raise ValueError(
                    "No structural data could be loaded. Is the input text in mol2 format?"
                )
            else:
                raise e

        # Infer residue PDB ID and name from substructure name
        mol2_df = cls._split_mol2_subst_names(mol2_df)

        # Format DataFrame
        mol2_df = cls._format_dataframe(mol2_df, verbose)

        return mol2_df

    @classmethod
    def _split_mol2_subst_names(cls, mol2_df):
        """
        Split "residue.subst_name" column values from mol2 file and add as
        "residue.name" and "residue.id" columns.

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
            lambda x: cls._split_mol2_subst_name(x["residue.subst_name"], x["atom.type"]),
            axis=1,
        )
        res_names = [res_name for res_name, res_id in result]
        res_ids = [res_id for res_name, res_id in result]

        mol2_df["residue.name"] = res_names
        mol2_df["residue.id"] = res_ids

        # Log a warning for all residues that cannot be converted to an integer
        for (res_id, subst_name), _ in mol2_df.groupby(
            ["residue.id", "residue.subst_name"], sort=False
        ):
            try:
                int(res_id)
            except ValueError:
                logger.warning(f"Suspicious residue ID: {res_id} (from {subst_name})")

        return mol2_df

    @classmethod
    def _split_mol2_subst_name(cls, subst_name, atom_type):
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
        - 1 ("residue.id").
        """

        # TODO in the future: Use regex!

        # Handle "residues" that are elements such as CA or MG.
        if subst_name[:2] == atom_type.upper():
            res_name = subst_name[:2]
            res_id = subst_name[2:]

        # These are amino acid, linkers, compounds, ...
        else:
            res_name = subst_name[:3]
            res_id = subst_name[3:]

        return res_name, res_id

    @classmethod
    def _format_dataframe(cls, dataframe, verbose=False):
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

                if column_dtype == "string":
                    dataframe.insert(len(dataframe.columns), column_name, None)
                elif column_dtype == "float32":
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
