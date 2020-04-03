"""
sequences.py

Utilities for sequence alignment
"""

import biotite
import biotite.sequence.align as align
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta

import warnings


def matrices(name: str) -> align.SubstitutionMatrix:

    """
    A SubstitutionMatrix maps each possible pairing of a symbol of a first alphabet with
    a symbol of a second alphabet to a score (int)

    Parameters
    ----------
    name: string
        Name of the matrix which is loaded from the internal matrix database.
        If the name of Substutition Matrix could not be found, the default SubstitutionMatrix
        will be BLOSUM62.

    Returns
    -------
    SubstitutionMatrix
        The class uses a 2-D (m x n) ndarray, where each element stores the score
        for a symbol pairing, indexed by the symbol codes of the respective symbols
        in an m-length alphabet 1 and an n-length alphabet 2

    """
    if name not in align.SubstitutionMatrix.list_db():
        raise ValueError(f"Substitution Matrix '{name}' could not be found.")
        matrix = align.SubstitutionMatrix.std_protein_matrix()
    else:
        alph = seq.ProteinSequence.alphabet
        matrix = align.SubstitutionMatrix(alph, alph, name)

    return matrix


def needleman_wunsch(seq1: str, seq2: str, matrix: str, gap: int) -> str:

    """
    Perform a global alignment, based on the Needleman-Wunsch algorithm

    Parameters
    ----------
    seq1,seq2: str
        The sequences to be aligned

    matrix: SubstitutionMatrix
        The substitution matrix used for scoring

    gap: int or (tuple, dtype=int)
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.

    Returns
    -------
    str
        An optimal alignment of two sequences
    """

    matrix = matrices(matrix)
    alignment = align.align_optimal(
        seq.ProteinSequence(seq1), seq.ProteinSequence(seq2), matrix, gap_penalty=gap
    )
    return alignment[0]


def smith_waterman(seq1: str, seq2: str, matrix: str, gap: int) -> str:

    """
    Perform a global alignment, based on the Smith-Waterman  algorithm

    Parameters
    ----------
    seq1,seq2: str
        The sequences to be aligned

    matrix: SubstitutionMatrix
        The substitution matrix used for scoring

    gap: int or (tuple, dtype=int)
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.

    Returns
    -------
    str
        A local alignment of two sequences
    """
    matrix = matrices(matrix)
    alignment = align.align_optimal(
        seq.ProteinSequence(seq1), seq.ProteinSequence(seq2), matrix, local=True, gap_penalty=gap
    )

    return alignment[0]


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
