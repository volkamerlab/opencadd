"""
opencadd.io.biopython

Defines classes that convert structural data into biopython structure objects.
"""

from pathlib import Path

from Bio.PDB.PDBParser import PDBParser


class PdbToBiopython:
    """
    Parse a pdb file into a biopython structure object.
    """

    @classmethod
    def from_file(cls, pdb_file):
        """
        Load structure as biopython structure from file.

        Parameters
        ----------
        pdb_file : pathlib.Path or str
            Path to pdb file.
        """

        pdb_file = Path(pdb_file)

        parser = PDBParser()
        structure = parser.get_structure("", pdb_file)

        return structure
