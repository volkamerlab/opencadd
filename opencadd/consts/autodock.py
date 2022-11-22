from typing import Literal, NamedTuple


class AutoDockAtomType(NamedTuple):
    """
    AutoDock-defined atom type with specified hydrogen-bonding properties.

    Attributes
    ----------
    name : str
        Name of the atom type according to AutoDock, e.g. H, HD, N, O, OA etc.
    description : str
        Short description of the atom type: element name, H-bonding role and type.
    h_bond_status : Literal[-1, 0, 1], Optional, default = 0
        +1 for H-bond acceptors (i.e. electron-pair donors),
        -1 for H-bond donors (i.e. electron-pair acceptors),
        0 for non H-bonding atom types.
    H_bond_count : int, Optional, default: 0
        Number of possible H-bonds for directionally H-bonding atoms, or -1 for spherically
        H-bonding atoms. For all non H-bonding atoms this must be set to 0 (default).

    See Also
    --------
    opencadd.consts.autodock.AutoDockAtomTypes
    """
    name: str
    description: str
    h_bond_status: Literal[-1, 0, 1] = 0
    h_bond_count: int = 0


class AutoDockAtomTypes(NamedTuple):
    """
    Collection of AutoDock-defined atom types and their properties,
    used in AutoDock software (e.g. AutoGrid4) and file formats (e.g. PDBQT).

    See Also
    --------
    opencadd.consts.autodock.AutoDockAtomType

    References
    ----------
    An official documentation could not be found; the list is taken from below link:
    https://mmb.irbbarcelona.org/gitlab/BioExcel/structureChecking/blob/5f07d82dc36d1f43733ae3b1ecd9f40aebe8b0a2/biobb_structure_checking/dat/autodock_atomtypes.dat
    """

    H: AutoDockAtomType = AutoDockAtomType("H", "Hydrogen, non H-bonding")
    H_D1: AutoDockAtomType = AutoDockAtomType("HD", "Hydrogen, H-bond donor, 1 bond", -1, 1)
    H_S: AutoDockAtomType = AutoDockAtomType("HS", "Hydrogen, H-bond donor, spherical", -1, -1)
    C_aliph: AutoDockAtomType = AutoDockAtomType("C", "Carbon, aliphatic, non H-bonding")
    C_arom: AutoDockAtomType = AutoDockAtomType("A", "Carbon, aromatic, non H-bonding")
    N: AutoDockAtomType = AutoDockAtomType("N", "Nitrogen, non H-bonding")
    NA: AutoDockAtomType = AutoDockAtomType("NA", "Nitrogen, H-bond acceptor, 1 bond", 1, 1)
    NS: AutoDockAtomType = AutoDockAtomType("NS", "Nitrogen, H-bond acceptor, spherical", 1, -1)
    OA: AutoDockAtomType = AutoDockAtomType("OA", "Oxygen, H-bond acceptor, 2 bonds", 1, 2)
    OS: AutoDockAtomType = AutoDockAtomType("OS", "Oxygen, H-bond acceptor, spherical", 1, -1)
    F: AutoDockAtomType = AutoDockAtomType("F", "Fluorine, non H-bonding")
    Mg: AutoDockAtomType = AutoDockAtomType("Mg", "Magnesium, non H-bonding")
    P: AutoDockAtomType = AutoDockAtomType("P", "Phosphorus, non H-bonding")
    SA: AutoDockAtomType = AutoDockAtomType("SA", "Sulfur, H-bond acceptor, 2 bonds", 1, 2)
    S: AutoDockAtomType = AutoDockAtomType("S", "Sulfur, non H-bonding")
    Cl: AutoDockAtomType = AutoDockAtomType("CL", "Chlorine, non H-bonding")
    Ca: AutoDockAtomType = AutoDockAtomType("CA", "Calcium, non H-bonding")
    Mn: AutoDockAtomType = AutoDockAtomType("Mn", "Manganese, non H-bonding")
    Fe: AutoDockAtomType = AutoDockAtomType("Fe", "Iron, non H-bonding")
    Zn: AutoDockAtomType = AutoDockAtomType("Zn", "Zinc, non H-bonding")
    Br: AutoDockAtomType = AutoDockAtomType("Br", "Bromine, non H-bonding")
    I: AutoDockAtomType = AutoDockAtomType("I", "Iodine, non H-bonding")
