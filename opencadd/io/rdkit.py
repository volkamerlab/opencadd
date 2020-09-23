"""
opencadd.io.rdkit

Defines classes that convert structural data into RDKit molecule objects.
"""

from rdkit import Chem
from rdkit.Chem import AllChem


class Mol2ToRdkitMol:
    """
    Parse a mol2 file or mol2 text into an RDKit molecule.
    """

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

        mol = Chem.MolFromMol2File(str(mol2_file))

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
