"""
Utilities for sequence alignment
"""
import logging

import numpy as np
import biotite.sequence.align as seq_align
import biotite.sequence as seq
import Bio

logger = logging.getLogger(__name__)


def matrices(name):

    """
    A SubstitutionMatrix maps each possible pairing of a symbol of a first alphabet with
    a symbol of a second alphabet to a score (int)

    Parameters
    ----------
    name: string
        Name of the matrix which is loaded from the internal matrix database.
        If the name of Substitution Matrix could not be found, the default SubstitutionMatrix
        will be BLOSUM62.

    Returns
    -------
    SubstitutionMatrix
        The class uses a 2-D (m x n) ndarray, where each element stores the score
        for a symbol pairing, indexed by the symbol codes of the respective symbols
        in an m-length alphabet 1 and an n-length alphabet 2

    """
    if name == "BLOSUM62":
        matrix = seq_align.SubstitutionMatrix.std_protein_matrix()
    else:
        alph = seq.ProteinSequence.alphabet
        matrix = seq_align.SubstitutionMatrix(alph, alph, name)

    return matrix


def sequence_alignment(seq1: str, seq2: str, matrix: str, gap: int, local: bool = False) -> str:

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

    local : bool, optional, default=False
        Whether to use local alignment (Smith-Waterman) or global (Needleman-Wunsch)

    Returns
    -------
    str
        An optimal alignment of two sequences
    """

    matrix = matrices(matrix)
    alignment = seq_align.align_optimal(
        seq.ProteinSequence(seq1),
        seq.ProteinSequence(seq2),
        matrix,
        gap_penalty=gap,
        local=local,
    )
    return alignment[0]


def fasta2select(
    fastafilename,
    ref_resids=None,
    ref_segids=None,
    target_resids=None,
    target_segids=None,
    backbone_selection="backbone or name CB",
):
    """Return selection strings that will select equivalent residues.

    The function takes two pre-aligned sequences in a FASTA file and
    constructs MDAnalysis selection strings of the common atoms. When
    these two strings are applied to the two different proteins they
    will generate AtomGroups of the aligned residues.


    Parameters
    ----------
    fastafilename : str, path to filename
        FASTA file with first sequence as reference and
        second the one to be aligned (ORDER IS IMPORTANT!)
    ref_resids : list (optional)
        sequence of resids as they appear in the reference structure
    target_resids : list (optional)
        sequence of resids as they appear in the target
    ref_segids : list (optional)
        sequence of segids as they appear in the reference structure
    target_segids : list (optional)
        sequence of segids as they appear in the target
    Returns
    -------
    dict
        dictionary with 'reference' and 'mobile' selection string

    """
    protein_gapped = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
    with open(fastafilename) as fasta:
        alignment = Bio.AlignIO.read(fasta, "fasta", alphabet=protein_gapped)

    nseq = len(alignment)
    if nseq != 2:
        raise ValueError("Only two sequences in the alignment can be processed.")

    orig_resids = [ref_resids, target_resids]
    orig_segids = [np.asarray(ref_segids), np.asarray(target_segids)]
    for iseq, a in enumerate(alignment):
        # need iseq index to change orig_resids
        if orig_resids[iseq] is None:
            # build default: assume consecutive numbering of all
            # residues in the alignment
            GAP = a.seq.alphabet.gap_char
            length = len(a.seq) - a.seq.count(GAP)
            orig_resids[iseq] = np.arange(1, length + 1)
        else:
            orig_resids[iseq] = np.asarray(orig_resids[iseq])
        # repeat for segment ids
        if orig_segids[iseq] is None:
            # build default: assume consecutive numbering of all
            # residues in the alignment
            GAP = a.seq.alphabet.gap_char
            length = len(a.seq) - a.seq.count(GAP)
            orig_segids[iseq] = np.full(length, "")
        else:
            orig_segids[iseq] = np.asarray(orig_segids[iseq])
    seq2resids = list(zip(orig_resids, orig_segids))

    def resid_factory(alignment, seq2resids):
        """Return a function that gives the resid for a position ipos in
        the nseq'th alignment.

        resid = resid_factory(alignment,seq2resids)
        r = resid(nseq,ipos)

        It is based on a look up table that translates position in the
        alignment to the residue number in the original
        sequence/structure.

        The first index of resid() is the alignmment number, the
        second the position in the alignment.

        seq2resids translates the residues in the sequence to resid
        numbers in the psf. In the simplest case this is a linear map
        but if whole parts such as loops are ommitted from the protein
        the seq2resids may have big gaps.

        Format: a tuple of two numpy arrays; the first array is for
        the reference, the second for the target, The index in each
        array gives the consecutive number of the amino acid in the
        sequence, the value the resid in the structure/psf.

        Note: assumes that alignments have same length and are padded if
        necessary.
        """
        # could maybe use Bio.PDB.StructureAlignment instead?
        nseq = len(alignment)
        t = np.zeros((nseq, alignment.get_alignment_length()), dtype=int)
        s = np.zeros((nseq, alignment.get_alignment_length()), dtype=object)
        for iseq, a in enumerate(alignment):
            GAP = a.seq.alphabet.gap_char
            indices = np.cumsum(np.where(np.array(list(a.seq)) == GAP, 0, 1)) - 1
            t[iseq, :] = seq2resids[iseq][0][indices]
            s[iseq, :] = seq2resids[iseq][1][indices]
            # -1 because seq2resid is index-1 based (resids start at 1)

        def resid(nseq, ipos, t=t, s=s):
            return t[nseq, ipos], s[nseq, ipos]

        return resid

    resid = resid_factory(alignment, seq2resids)
    res_list = []  # collect individual selection string

    # should be the same for both seqs
    GAP = alignment[0].seq.alphabet.gap_char
    if GAP != alignment[1].seq.alphabet.gap_char:
        raise ValueError("Different gap characters in sequence 'target' and 'mobile'.")
    for ipos in range(alignment.get_alignment_length()):
        aligned = list(alignment[:, ipos])
        if GAP in aligned:
            continue  # skip residue if it's not matching any
        template = "( resid {:d} and segid {:s}" + f" and ( {backbone_selection} ) )"
        res_list.append([template.format(*resid(iseq, ipos)) for iseq in range(nseq)])

    sel = np.array(res_list).transpose()

    ref_selection = " or ".join(sel[0])
    target_selection = " or ".join(sel[1])
    return {"reference": ref_selection, "mobile": target_selection}
