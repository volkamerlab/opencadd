"""
opencadd.io.api

Defines the opencadd.oi API.
"""

from pathlib import Path

from . import dataframes, rdkit


class DataFrame:
    """
    Core object to load structures as DataFrame.
    """

    @classmethod
    def from_file(cls, filepath, verbose=False):
        """
        Load structures as DataFrame from file.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to structure file: pdb and mol2 files only.
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).
        """

        filepath = Path(filepath)
        format = filepath.suffix[1:]

        if format == "mol2":
            mol2_to_df = dataframes.Mol2ToDataFrame()
            dataframe = mol2_to_df.from_file(filepath, verbose)
        elif format == "pdb":
            pdb_to_df = dataframes.PdbToDataFrame()
            dataframe = pdb_to_df.from_file(filepath, verbose)
        else:
            raise ValueError(f"The {format} format is not suppored or invalid.")

        return dataframe

    @classmethod
    def from_text(cls, text, format, verbose=False):
        """
        Load structures as DataFrame from text (file content as string).

        Parameters
        ----------
        text : str
            Structure file content as string.
        format : str
            Structure format: "pdb" or "mol2".
        verbose : bool
            Show only default columns (False) or additionally input-format specific columns (True).
        """

        if format == "mol2":
            mol2_to_df = dataframes.Mol2ToDataFrame()
            dataframe = mol2_to_df.from_text(text, verbose)
        elif format == "pdb":
            pdb_to_df = dataframes.PdbToDataFrame()
            dataframe = pdb_to_df.from_text(text, verbose)
        else:
            raise ValueError(f"The {format} format is not suppored or invalid.")

        return dataframe


class RdkitMol:
    """
    Core object to load structures as RDKit molecule.
    """

    @classmethod
    def from_file(cls, filepath, compute2d=True):
        """
        Load structures as DataFrame from file.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to structure file: pdb and mol2 files only.
        compute2d : bool
            Compute 2D coordinates for ligand (default).
        """

        filepath = Path(filepath)
        format = filepath.suffix[1:]

        if format == "mol2":
            mol2_to_rmol = rdkit.Mol2ToRdkitMol()
            rmol = mol2_to_rmol.from_file(filepath, compute2d)
        else:
            raise ValueError(f"The {format} format is not suppored or invalid.")

        return rmol

    @classmethod
    def from_text(cls, text, format, compute2d=True):
        """
        Load structures as DataFrame from text (file content as string).

        Parameters
        ----------
        text : str
            Structure file content as string.
        format : str
            Structure format: "pdb" or "mol2".
        compute2d : bool
            Compute 2D coordinates for ligand (default).
        """

        if format == "mol2":
            mol2_to_rmol = rdkit.Mol2ToRdkitMol()
            rmol = mol2_to_rmol.from_text(text, compute2d)
        else:
            raise ValueError(f"The {format} format is not suppored or invalid.")

        return rmol