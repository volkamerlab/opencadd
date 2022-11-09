"""
Functions and routines for easier communication with PyMOL.
"""


# Standard library
from pathlib import Path
# 3rd party
from pymol import cmd
import numpy as np


# scans each grid point for pseudoatoms inside HBond distance (input parameter)
def computeHBondEnvironment(
        coordinates: np.ndarray,
        hbond_len: float,
        filepath_pdb: Path,
):
    """
    Calculate # HBond donors and acceptors per grid points within %f A" % (hbond_len)

    Parameters
    ----------
    coordinates
    hbond_len
    filepath_pdb

    Returns
    -------

    """

    # variable names flipped, donor-acceptor ist Ansichtssache je nach von Protein oder Ligand
    # betrachtet
    # First element is list of acceptor atoms, and second element is list of donors
    acceptor_donor_atoms = {
        "ALA": [["N"], ["O"]],
        "ARG": [["NE", "NH1", "NH2", "N"], ["O"]],
        "ASN": [["NH2", "N"], ["OD1", "O"]],
        "ASP": [["OD2", "N"], ["OD1", "OD2", "O"]],
        "CYS": [["N"], ["O"]],
        "GLU": [["OE1", "N"], ["OE1", "OE2", "O"]],
        "GLN": [["NE2", "N"], ["OE1", "O"]],
        "GLY": [["N"], ["O"]],
        "HIS": [["ND1", "N"], ["NE2", "O"]],
        "ILE": [["N"], ["O"]],
        "LEU": [["N"], ["O"]],
        "LYS": [["NZ", "N"], ["O"]],
        "MET": [["N"], ["O"]],
        "PHE": [["N"], ["O"]],
        "PRO": [["N"], ["O"]],
        "SER": [["OG", "N"], ["OG", "O"]],
        "THR": [["OG1", "N"], ["OG1", "O"]],
        "TRP": [["NE1", "N"], ["O"]],
        "TYR": [["OH", "N"], ["OH", "O"]],
        "VAL": [["N"], ["O"]],
        "HOH": [["O"], ["H"]],
        "ZN": [["ZN"], ["ZN"]],
    }

    cmd.reinitialize()
    cmd.load(filepath_pdb, object="proteinfile")
    donors_count = np.zeros(coordinates.shape[0], dtype=np.uint8)
    acceptors_count = np.zeros(coordinates.shape[0], dtype=np.uint8)
    for idx, coord in enumerate(coordinates):
        cmd.pseudoatom("tmpPoint", pos=coord)
        # radius of circle for energy density
        cmd.select("proteinAtomshydrophil", f"(tmpPoint around {hbond_len}) and proteinfile")
        atoms_data = {"stored_names": []}
        cmd.iterate_state(
            1, "proteinAtomshydrophil", "stored_names.append([resn, name])", space=atoms_data
        )
        cmd.delete("tmpPoint")
        cmd.delete("proteinAtomshydrophil")
        atoms_data = atoms_data["stored_names"]
        # calculates number of donor/acceptor atoms in neighbourhood
        for res_name, atom_type in atoms_data:
            if atom_type in acceptor_donor_atoms[res_name][0]:
                acceptors_count[idx] += 1
            if atom_type in acceptor_donor_atoms[res_name][1]:
                donors_count[idx] += 1
    return acceptors_count, donors_count
