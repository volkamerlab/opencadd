import biotite
import biotite.database.rcsb as rcsb
import biotite.structure.io.pdb as pdb
import biotite.sequence as seq
import biotite.structure.io as strucio
import biotite.structure as struc

from needleman_wunsch import needleman_wunsch
from smith_waterman import smith_waterman


def biotite_aminoSeq(name: str) -> str:

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

    pdb_file_path = rcsb.fetch(name[0 : len(name) - 4], "pdb", biotite.temp_dir())
    file = pdb.PDBFile()
    file.read(pdb_file_path)
    atomarraystack = strucio.load_structure(pdb_file_path)
    ca_atoms = atomarraystack[atomarraystack.atom_name == "CA"]
    residuen = ca_atoms.res_name
    sequence = list(map(seq.ProteinSequence.convert_letter_3to1, residuen))

    return "".join([str(elem) for elem in sequence])


# print(type(residuen).__name__)


def get_alignment(ref: str, mobile: str, methods: bool, matrix: str, gap: int) -> str:

    """
    Perform a global alignment, based on the the Needleman-Wunsch algorithm
    or a local alignment, based on the Smith-Waterman algorithm

    Parameters
    ----------
    ref,mobile: string
        string of sequences

    methods: bool
        true for local alignment, otherwise a global alignment is performed

    matrix: string
        The substitution matrix used for scoring

    gap: int or (tuple, dtype=int)
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.

    Returns
    -------
    string
    An alignment of two sequences

    """

    seq1, seq2 = biotite_aminoSeq(ref), biotite_aminoSeq(mobile)

    if methods == True:
        alignment = smith_waterman(seq1, seq2, matrix, gap)
    else:
        alignment = needleman_wunsch(seq1, seq2, matrix, gap)

    return alignment


# print(get_alignment("4u3y.pdb", "4u40.pdb", False, "PAM250", (-5, -15)))
