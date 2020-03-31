import biotite
import biotite.database.rcsb as rcsb
import biotite.structure.io.pdb as pdb
import biotite.sequence as seq
import biotite.structure.io as strucio
import biotite.structure as struc
import MDAnalysis as mda
import biotite.sequence.io.fasta as fasta

from structuralalignment.alignment.needleman_wunsch import needleman_wunsch
from structuralalignment.alignment.smith_waterman import smith_waterman


def biotite_amino_seq(name: str) -> str:

    """
    1. Download structure file from RCSB
    2. Open the file and get the amino acid sequence

    Parameter
    ---------
    name: string
          name of the file in string

    Return
    ------
    string
        one-letter amino acid sequence of the file

    """
    universe = mda.Universe(name)
    ca_atoms = universe.select_atoms("name CA")
    sequence = ca_atoms.residues.sequence()
    return str(sequence.seq)


def get_alignment_fasta(fasta_file):

    """
    Get an alignment from a FastaFile instance.

    Parameter
    ---------
    fasta_file: FastaFile

    Returns
    -------
    alignment: Alignment
    """
    sequence = fasta.get_alignment(fasta_file)
    return sequence


def get_alignment(
    ref: str, mobile: str, local: bool = False, matrix: str = "BLOSUM62", gap: int = -10
) -> str:

    """
    Perform a global alignment, based on the the Needleman-Wunsch algorithm
    or a local alignment, based on the Smith-Waterman algorithm

    Parameters
    ----------
    ref,mobile: string
        string of sequences

    methods: bool, optional
        true for local alignment, otherwise a global alignment is performed
        (Default: False)

    matrix: string, optional
        The substitution matrix used for scoring
        (Default: BLOSUM62)

    gap: int or (tuple, dtype=int), optional
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.
        (Default: -10)

    Returns
    -------
    string
    An alignment of two sequences

    Examples
    --------
    >>> get_alignment("xxxx.pdb", "xxxx.pdb", False, "PAM250", (-5, -15))

    RKKSLVDIDLSSLRDP
    R-K-I-DLS-S-LRDP

    """

    seq1, seq2 = biotite_amino_seq(ref), biotite_amino_seq(mobile)

    if local is True:
        alignment = smith_waterman(seq1, seq2, matrix, gap)
    else:
        alignment = needleman_wunsch(seq1, seq2, matrix, gap)

    return alignment


if __name__ == "__main__":
    print(get_alignment_fasta("4u40.fasta.txt"))
# print(get_alignment("4u3y.pdb", "4u40.pdb"))

