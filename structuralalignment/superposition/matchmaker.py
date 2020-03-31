import biotite
import biotite.structure as struc
from structuralalignment.alignment.base import get_alignment
from structuralalignment.alignment.base import biotite_amino_seq

import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align


def get_trace_atoms(ref: str, mobile: str, alignment: str) -> mda.AtomGroup:

    """
    Get the trace of the aligned residues in the sequences and use it
    to map corresponding residues in their structures

    Parameters
    ----------
    ref,mobile: string
        name of file

    alignment: string
        an alignment of two sequences

    Returns
    -------
    AtomGroup with atoms
        AtomGroup with information about the atoms such as names,
        indices, or the coordinates in the positions attribute

    Examples
    --------
    >>> a = get_alignment("xxxx.pdb", "xxxx.pdb", True, "PAM250", (-5, -15))
    get_trace_atoms("xxxx.pdb", "xxxx.pdb", a)

    [<Atom 1: N of type N of resname SER, resid 17 and segid A and altLoc>,....]

    """
    ref_universe = mda.Universe(ref)
    mobile_universe = mda.Universe(mobile)
    trace = alignment.trace
    trace = trace[~(trace == -1).any(axis=1)]
    aref = ref_universe.residues[trace[:, 0]]
    amob = mobile_universe.residues[trace[:, 1]]

    return aref.atoms, amob.atoms


def align_protein(ref: mda.AtomGroup, mobile: mda.AtomGroup) -> (float, float):

    """
    Perform a spatial superposition by minimizing the RMSD.
    Spatially align the group of atoms mobile to reference by doing a RMSD fit on select atoms.

    Parameters
    ----------
    ref,mobile: AtomGroup with Atoms

    Returns
    -------
    old_rmsd (float) – RMSD before spatial alignment
    new_rmsd (float) – RMSD after spatial alignment

    """
    return mda_align.alignto(ref, mobile, strict=False, select="protein and name CA")


def create_pdbfile(ref: str, mobile: str, alignment: str) -> str:

    """
    Create a new PDB-file

    Parameters
    ----------
    ref,mobile: string
        name of file

    alignment: string
        an alignment of two sequences

    Returns
    -------

    PDF-file
        a new PDF-file

    string
        The PDB-file was successfully created.

    """

    aref, amob = get_trace_atoms(ref, mobile, alignment)
    mda_align.alignto(aref, amob, strict=False, select="protein and name CA")
    with mda.Writer("protein.pdb", multiframe=True) as pdb:
        pdb.write(aref)
        pdb.write(amob)
    return "The PDB-file was successfully created. "


if __name__ == "__main__":
    a = get_alignment("4u3y.pdb", "4u40.pdb", True, "PAM250", (-5, -15))
    b1, b2 = get_trace_atoms("4u3y.pdb", "4u40.pdb", a)
    # print(align_protein(b1, b2))
    print(create_pdbfile("4u3y.pdb", "4u40.pdb", a))
