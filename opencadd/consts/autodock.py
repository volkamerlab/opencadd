from typing import Literal, NamedTuple


class AtomType(NamedTuple):
    """
    AutoDock-defined atom type with specified hydrogen-bonding properties.

    Attributes
    ----------
    name : str
        Name of the atom type according to AutoDock, e.g. H, HD, N, O, OA etc.
    description : str
        Short description of the atom type: element name, H-bonding role and type.
    hbond_status : Literal[-1, 0, 1], Optional, default = 0
        +1 for H-bond acceptors (i.e. electron-pair donors),
        -1 for H-bond donors (i.e. electron-pair acceptors),
        0 for non H-bonding atom types.
    hbond_count : int, Optional, default: 0
        Number of possible H-bonds for directionally H-bonding atoms, or -1 for spherically
        H-bonding atoms. For all non H-bonding atoms this must be set to 0 (default).

    See Also
    --------
    opencadd.consts.autodock.AutoDockAtomTypes
    """
    name: str
    description: str
    hbond_status: Literal[-1, 0, 1] = 0
    hbond_count: int = 0


class AtomTypes(NamedTuple):
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

    H: AtomType = AtomType("H", "Hydrogen, non H-bonding")
    H_D1: AtomType = AtomType("HD", "Hydrogen, H-bond donor, 1 bond", -1, 1)
    H_S: AtomType = AtomType("HS", "Hydrogen, H-bond donor, spherical", -1, -1)
    C_aliph: AtomType = AtomType("C", "Carbon, aliphatic, non H-bonding")
    C_arom: AtomType = AtomType("A", "Carbon, aromatic, non H-bonding")
    N: AtomType = AtomType("N", "Nitrogen, non H-bonding")
    NA: AtomType = AtomType("NA", "Nitrogen, H-bond acceptor, 1 bond", 1, 1)
    NS: AtomType = AtomType("NS", "Nitrogen, H-bond acceptor, spherical", 1, -1)
    OA: AtomType = AtomType("OA", "Oxygen, H-bond acceptor, 2 bonds", 1, 2)
    OS: AtomType = AtomType("OS", "Oxygen, H-bond acceptor, spherical", 1, -1)
    F: AtomType = AtomType("F", "Fluorine, non H-bonding")
    Mg: AtomType = AtomType("Mg", "Magnesium, non H-bonding")
    P: AtomType = AtomType("P", "Phosphorus, non H-bonding")
    SA: AtomType = AtomType("SA", "Sulfur, H-bond acceptor, 2 bonds", 1, 2)
    S: AtomType = AtomType("S", "Sulfur, non H-bonding")
    Cl: AtomType = AtomType("CL", "Chlorine, non H-bonding")
    Ca: AtomType = AtomType("CA", "Calcium, non H-bonding")
    Mn: AtomType = AtomType("Mn", "Manganese, non H-bonding")
    Fe: AtomType = AtomType("Fe", "Iron, non H-bonding")
    Zn: AtomType = AtomType("Zn", "Zinc, non H-bonding")
    Br: AtomType = AtomType("Br", "Bromine, non H-bonding")
    I: AtomType = AtomType("I", "Iodine, non H-bonding")
