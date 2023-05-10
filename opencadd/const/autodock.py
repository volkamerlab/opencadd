from typing import Literal
from enum import Enum


class AtomType(Enum):
    """
    Collection of AutoDock-defined atom types and their properties,
    used in AutoDock software (e.g. AutoGrid4) and file formats (e.g. PDBQT).

    References
    ----------
    An official documentation could not be found; the list is taken from below link:
    https://mmb.irbbarcelona.org/gitlab/BioExcel/structureChecking/blob/5f07d82dc36d1f43733ae3b1ecd9f40aebe8b0a2/biobb_structure_checking/dat/autodock_atomtypes.dat
    """
    H = ("Hydrogen, non H-bonding", 0, 0)
    HD = ("Hydrogen, H-bond donor, 1 bond", -1, 1)
    HS = ("Hydrogen, H-bond donor, spherical", -1, -1)
    C = ("Carbon, aliphatic, non H-bonding", 0, 0)
    A = ("Carbon, aromatic, non H-bonding", 0, 0)
    N = ("Nitrogen, non H-bonding", 0, 0)
    NA = ("Nitrogen, H-bond acceptor, 1 bond", 1, 1)
    NS = ("Nitrogen, H-bond acceptor, spherical", 1, -1)
    OA = ("Oxygen, H-bond acceptor, 2 bonds", 1, 2)
    OS = ("Oxygen, H-bond acceptor, spherical", 1, -1)
    F = ("Fluorine, non H-bonding", 0, 0)
    Mg = ("Magnesium, non H-bonding", 0, 0)
    P = ("Phosphorus, non H-bonding", 0, 0)
    SA = ("Sulfur, H-bond acceptor, 2 bonds", 1, 2)
    S = ("Sulfur, non H-bonding", 0, 0)
    Cl = ("Chlorine, non H-bonding", 0, 0)
    Ca = ("Calcium, non H-bonding", 0, 0)
    Mn = ("Manganese, non H-bonding", 0, 0)
    Fe = ("Iron, non H-bonding", 0, 0)
    Zn = ("Zinc, non H-bonding", 0, 0)
    Br = ("Bromine, non H-bonding", 0, 0)
    I = ("Iodine, non H-bonding", 0, 0)

    @property
    def description(self) -> str:
        """
        Short description of the atom type: element name, H-bonding role and count.
        """
        return self.value[0]

    @property
    def hbond_status(self) -> Literal[-1, 0, 1]:
        """
        Whether the atom type is an H-bond donor, H-bond acceptor, or non H-bonding.
        Returns
        -------
        Literal[-1, 0, 1]
            +1 for H-bond acceptors (i.e. electron-pair donors),
            -1 for H-bond donors (i.e. electron-pair acceptors),
            0 for non H-bonding atom types.
        """
        return self.value[1]

    @property
    def hbond_count(self) -> int:
        """
        Number of possible H-bonds for directionally H-bonding atoms, or -1 for spherically
        H-bonding atoms, and 0 for non H-bonding atoms.
        """
        return self.value[2]
