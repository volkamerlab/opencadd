"""
opencadd.io.rdkit

Defines classes that convert structural data into RDKit molecule objects.
"""

import logging

from rdkit import Chem
from rdkit.Chem import AllChem

from .core import _Base

logger = logging.getLogger(__name__)


class Rdkit(_Base):
    """
    Parse a structure as an RDKit molecule.
    """

    @classmethod
    def from_file(cls, filepath, compute2d=True):
        """
        Load structures as RDKit molecule from file.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to structure file: pdb and mol2 files only.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Structure as RDKit molecule object.
        """

        filepath = cls._convert_filepath(filepath)

        if filepath.suffix == ".mol2":
            return cls._from_mol2_file(filepath, compute2d)
        elif filepath.suffix == ".pdb":
            return cls._from_pdb_file(filepath, compute2d)
        else:
            raise ValueError(f"The {filepath.suffix} format is not supported or invalid.")

    @classmethod
    def from_text(cls, text, ext, compute2d=True):
        """
        Load structures as RDKit molecule from text (file content as string).

        Parameters
        ----------
        text : str
            Structure file content as string.
        ext : str
            Structure format: "pdb" or "mol2".
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Structure as RDKit molecule object.
        """

        if ext == "mol2":
            return cls._from_mol2_text(text, compute2d)
        else:
            raise ValueError(f"The {ext} format is not supported or invalid.")

    @classmethod
    def _from_mol2_file(cls, mol2_file, compute2d=True):
        """
        Get structural data from mol2 file as RDKit molecule.

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

    @classmethod
    def _from_mol2_text(cls, mol2_text, compute2d=True):
        """
        Get structural data from mol2 text as RDKit molecule.

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

    @classmethod
    def _from_pdb_file(cls, pdb_file, compute2d=True):
        """
        Get structural data from pdb file as RDKit molecule.

        Parameters
        ----------
        pdb_file : pathlib.Path or str
            Path to pdb file.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Molecule.
        """

        mol = Chem.MolFromPDBFile(str(pdb_file))

        if compute2d:
            AllChem.Compute2DCoords(mol)

        return mol

    @classmethod
    def _from_pdb_text(cls, pdb_text, compute2d=True):
        """
        Get structural data from pdb text as RDKit molecule.

        Parameters
        ----------
        mol2_text : str
            Pdb file content from KLIFS database.
        compute2d : bool
            Compute 2D coordinates for ligand (default).

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            Molecule.
        """

        mol = Chem.MolFromPDBBlock(pdb_text)

        if compute2d:
            AllChem.Compute2DCoords(mol)

        return mol
