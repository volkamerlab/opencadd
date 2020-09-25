"""
opencadd.io.biopython

Defines classes that convert structural data into biopython structure objects.
"""

import logging

from Bio.PDB.PDBParser import PDBParser

from .core import _Base

logger = logging.getLogger(__name__)


class Biopython(_Base):
    """
    Parse a structure as a biopython structure object.
    """

    @classmethod
    def from_file(cls, filepath):
        """
        Load structures as biopython structure object from file.

        Parameters
        ----------
        filepath : str or pathlib.Path
            Path to structure file: pdb files only.

        Returns
        -------
        Bio.PDB.Structure.Structure
            Structure as biopython structure object.
        """

        filepath = cls._convert_filepath(filepath)

        if filepath.suffix == ".pdb":
            return cls._from_pdb_file(filepath)
        else:
            raise ValueError(f"The {format} format is not supported or invalid.")

    @classmethod
    def _from_pdb_file(cls, pdb_file):
        """
        Load structure as biopython structure from file.

        Parameters
        ----------
        pdb_file : pathlib.Path
            Path to pdb file.
        """

        parser = PDBParser()
        structure = parser.get_structure("", pdb_file)

        return structure
