"""
Core objects used across the library
"""

from MDAnalysis import Universe


class Structure(Universe):
    """
    Core object to load structures with.

    Thin wrapper around MDAnalysis.Universe objects
    """

    @classmethod
    def from_pdbid(cls, pdbid):
        import mmtf

        return cls(mmtf.fetch(pdbid))

    @classmethod
    def from_string(cls, pdbid_or_path):
        import os

        if os.path.isfile(pdbid_or_path):
            return cls(pdbid_or_path)
        return cls.from_pdbid(pdbid_or_path)

    def write(self, *args, **kwargs):
        return self.atoms.write(*args, **kwargs)
