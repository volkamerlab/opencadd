"""
Utilities for sequence alignment
"""
import logging
import os

import numpy as np
import biotite.sequence.align as align
import biotite.sequence as seq
import biotite.sequence.io.fasta as fasta
import Bio

logger = logging.getLogger(__name__)


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
    if name == "BLOSUM62":
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
        a local alignment of two sequences
    """
    matrix = matrices(matrix)
    alignment = align.align_optimal(
        seq.ProteinSequence(seq1), seq.ProteinSequence(seq2), matrix, local=True, gap_penalty=gap,
    )

    return alignment[0]


def get_alignment_fasta(fasta_file):
    """
    Get an alignment from a FastaFile instance.

    Parameters
    ----------
    fasta_file: FastaFile

    Returns
    -------
    alignment: Alignment
    """
    sequence = fasta.get_alignment(fasta_file)
    return sequence


def fasta2select(
    fastafilename,
    is_aligned=False,
    ref_resids=None,
    ref_segids=None,
    target_resids=None,
    target_segids=None,
    ref_offset=0,
    target_offset=0,
    verbosity=3,
    alnfilename=None,
    treefilename=None,
    backbone_selection="backbone or name CB",
    clustalw="clustalw2",
):
    """Return selection strings that will select equivalent residues.

    The function aligns two sequences provided in a FASTA file and
    constructs MDAnalysis selection strings of the common atoms. When
    these two strings are applied to the two different proteins they
    will generate AtomGroups of the aligned residues.

    `fastafilename` contains the two un-aligned sequences in FASTA
    format. The reference is assumed to be the first sequence, the
    target the second. ClustalW_ produces a pairwise
    alignment (which is written to a file with suffix ``.aln``).  The
    output contains atom selection strings that select the same atoms
    in the two structures.

    Unless `ref_offset` and/or `target_offset` are specified, the resids
    in the structure are assumed to correspond to the positions in the
    un-aligned sequence, namely the first residue has resid == 1.

    In more complicated cases (e.g., when the resid numbering in the
    input structure has gaps due to missing parts), simply provide the
    sequence of resids as they appear in the topology in `ref_resids` or
    `target_resids`, e.g. ::

       target_resids = [a.resid for a in trj.select_atoms('name CA')]

    (This translation table *is* combined with any value for
    `ref_offset` or `target_offset`!)

    Parameters
    ----------
    fastafilename : str, path to filename
        FASTA file with first sequence as reference and
        second the one to be aligned (ORDER IS IMPORTANT!)
    is_aligned : bool (optional)
        ``False`` (default)
            run clustalw for sequence alignment;
        ``True``
            use the alignment in the file (e.g. from STAMP) [``False``]
    ref_offset : int (optional)
        add this number to the column number in the FASTA file
        to get the original residue number, default: 0
    target_offset : int (optional)
        add this number to the column number in the FASTA file
        to get the original residue number, default: 0
    ref_resids : str (optional)
        sequence of resids as they appear in the reference structure
    target_resids : str (optional)
        sequence of resids as they appear in the target
    alnfilename : str (optional)
        filename of ClustalW alignment (clustal format) that is
        produced by *clustalw* when *is_aligned* = ``False``.
        default ``None`` uses the name and path of *fastafilename* and
        subsititutes the suffix with '.aln'.
    treefilename: str (optional)
        filename of ClustalW guide tree (Newick format);
        if default ``None``  the the filename is generated from *alnfilename*
        with the suffix '.dnd' instead of '.aln'
    clustalw : str (optional)
        path to the ClustalW (or ClustalW2) binary; only
        needed for `is_aligned` = ``False``, default: "ClustalW2"

    Returns
    -------
    select_dict : dict
        dictionary with 'reference' and 'mobile' selection string
        that can be used immediately in :class:`AlignTraj` as
        ``select=select_dict``.


    See Also
    --------
    :func:`sequence_alignment`, which does not require external
    programs.


    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/

    """
    protein_gapped = Bio.Alphabet.Gapped(Bio.Alphabet.IUPAC.protein)
    if is_aligned:
        logger.debug("Using provided alignment %s", fastafilename)
        with open(fastafilename) as fasta:
            alignment = Bio.AlignIO.read(fasta, "fasta", alphabet=protein_gapped)
    else:
        if alnfilename is None:
            filepath, ext = os.path.splitext(fastafilename)
            alnfilename = filepath + ".aln"
        if treefilename is None:
            filepath, ext = os.path.splitext(alnfilename)
            treefilename = filepath + ".dnd"
        run_clustalw = Bio.Align.Applications.ClustalwCommandline(
            clustalw,
            infile=fastafilename,
            type="protein",
            align=True,
            outfile=alnfilename,
            newtree=treefilename,
        )
        logger.debug("Aligning sequences in %(fastafilename)r with %(clustalw)r.", vars())
        logger.debug("ClustalW commandline: %r", str(run_clustalw))
        try:
            stdout, stderr = run_clustalw()
        except:
            logger.exception("ClustalW %(clustalw)r failed", vars())
            logger.info("(You can get clustalw2 from http://www.clustal.org/clustal2/)")
            raise
        with open(alnfilename) as aln:
            alignment = Bio.AlignIO.read(aln, "clustal", alphabet=protein_gapped)
        logger.debug("Using clustalw sequence alignment {0!r}".format(alnfilename))
        logger.debug("ClustalW Newick guide tree was also produced: {0!r}".format(treefilename))

    nseq = len(alignment)
    if nseq != 2:
        raise ValueError("Only two sequences in the alignment can be processed.")

    # implict assertion that we only have two sequences in the alignment
    orig_resids = [ref_resids, target_resids]
    offsets = [ref_offset, target_offset]
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
    # add offsets to the sequence <--> resid translation table
    seq2resids = [
        (resids + offset, segids)
        for resids, offset, segids in zip(orig_resids, offsets, orig_segids)
    ]
    del orig_resids
    del offsets

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
    # could collect just resid and type (with/without CB) and
    # then post-process and use ranges for continuous stretches, eg
    # ( resid 1:35 and ( backbone or name CB ) ) or ( resid 36 and backbone )

    # should be the same for both seqs
    GAP = alignment[0].seq.alphabet.gap_char
    if GAP != alignment[1].seq.alphabet.gap_char:
        raise ValueError("Different gap characters in sequence 'target' and 'mobile'.")
    for ipos in range(alignment.get_alignment_length()):
        aligned = list(alignment[:, ipos])
        if GAP in aligned:
            continue  # skip residue
        template = "( resid {:d} and segid {:s}" + f" and ( {backbone_selection} ) )"

        res_list.append([template.format(*resid(iseq, ipos)) for iseq in range(nseq)])

    sel = np.array(res_list).transpose()

    ref_selection = " or ".join(sel[0])
    target_selection = " or ".join(sel[1])
    return {"reference": ref_selection, "mobile": target_selection}
