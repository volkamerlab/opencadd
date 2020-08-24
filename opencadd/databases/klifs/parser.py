"""
parser.py

Defines classes to parse data with different input/output.
"""

from pathlib import Path


from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

MOL2_COLUMNS = {
    "n_cols_10": {
        0: ("atom.id", int),
        1: ("atom.name", str),
        2: ("atom.x", float),
        3: ("atom.y", float),
        4: ("atom.z", float),
        5: ("atom.type", str),
        6: ("residue.subst_id", int),
        7: ("residue.subst_name", str),
        8: ("atom.charge", float),
        9: ("atom.backbone", str),
    },
    "n_cols_9": {
        0: ("atom.id", int),
        1: ("atom.name", str),
        2: ("atom.x", float),
        3: ("atom.y", float),
        4: ("atom.z", float),
        5: ("atom.type", str),
        6: ("residue.subst_id", int),
        7: ("residue.subst_name", str),
        8: ("atom.charge", float),
    },
}


class PdbToDataFrame:
    """
    Bla TODO
    """

    def __init__(self):
        pass

    @staticmethod
    def from_text(pdb_text):
        """
        Get structural data from pdb text.

        Parameters
        ----------
        pdb_text : str
            Pdb file content from KLIFS database.

        Returns
        -------
        dict of pandas.DataFrame
            Structural data
        """

        ppdb = PandasPdb()
        pdb_dict = ppdb._construct_df(pdb_text.splitlines(True))

        # Streamline output with mol2 output? TODO

        return pdb_dict


class Mol2ToDataFrame:
    """
    Bla TODO
    """

    def __init__(self):

        pass

    def from_file(self, mol2_file):
        """
        Get structural data from mol2 file.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
        Path to mol2 file.

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

        mol2_file = Path(mol2_file)
        if not mol2_file.exists():
            raise FileNotFoundError(f"File {mol2_file} does not exist.")

        pmol = PandasMol2()

        try:
            mol2_df = pmol.read_mol2(str(mol2_file), columns=MOL2_COLUMNS["n_cols_10"]).df

        except ValueError:
            mol2_df = pmol.read_mol2(str(mol2_file), columns=MOL2_COLUMNS["n_cols_9"]).df
            # Add column to return same column names in both cases (i.e. try-except)
            mol2_df.insert(9, "atom.backbone", np.nan)

        return self._split_mol2_subst_names(mol2_df)

    def from_text(self, mol2_text):
        """
        Get structural data from mol2 text.

        Parameters
        ----------
        mol2_text : str
            Mol2 file content from KLIFS database.

        Returns
        -------
        pandas.DataFrame
            Structural data.
        """

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
            # Add column to return same column names in both cases (i.e. try-except)
            mol2_df.insert(9, "atom.backbone", np.nan)

        return self._split_mol2_subst_names(mol2_df)

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


class Mol2ToRdkitMol:
    """
    Bla TODO
    """

    def __init__(self):
        pass

    @staticmethod
    def from_file(mol2_file, compute2d=True):
        """
        Get structural data from mol2 file.

        Parameters
        ----------
        mol2_file : pathlib.Path or str
        Path to mol2 file.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Molecule.
        """

        mol = Chem.MolFromMol2File(mol2_file)

        if compute2d:
            AllChem.Compute2DCoords(mol)

        return mol

    @staticmethod
    def from_text(mol2_text, compute2d=True):
        """
        Get structural data from mol2 text.

        Parameters
        ----------
        mol2_text : str
            Mol2 file content from KLIFS database.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Molecule.
        """

        mol = Chem.MolFromMol2Block(mol2_text)

        if compute2d:
            AllChem.Compute2DCoords(mol)

        return mol
