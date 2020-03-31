import biotite
import biotite.sequence as seq
import biotite.sequence.align as align
from structuralalignment.alignment.matrices import matrices


def needleman_wunsch(seq1: str, seq2: str, matrix: str, gap: int) -> str:

    """
    Perform a global alignment, based on the Needleman-Wunsch algorithm

    Parameters
    ----------
    seq1,seq2: string
        The sequences to be aligned

    matrix: SubstitutionMatrix
        The substitution matrix used for scoring

    gap: int or (tuple, dtype=int)
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.

    Return
    ------
    string
    An optimal alignment of two sequences
    """

    matrix = matrices(matrix)
    alignment = align.align_optimal(
        seq.ProteinSequence(seq1), seq.ProteinSequence(seq2), matrix, gap_penalty=gap,
    )
    return alignment[0]

