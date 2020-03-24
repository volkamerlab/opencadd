import biotite
import biotite.sequence.align as align
import biotite.sequence as seq
import warnings


def matrices(name: str) -> align.SubstitutionMatrix:

    """
    A SubstitutionMatrix maps each possible pairing of a symbol of a first alphabet with
    a symbol of a second alphabet to a score (int)

    Parameter
    ---------
    name: string
        Name of the matrix which is loaded from the internal matrix database.
        If the name of Substutition Matrix could not be found, the default SubstitutionMatrix
        will be BLOSUM62.

    Return
    ------
    SubstitutionMatrix
        The class uses a 2-D (m x n) ndarray, where each element stores the score
        for a symbol pairing, indexed by the symbol codes of the respective symbols
        in an m-length alphabet 1 and an n-length alphabet 2

    """
    if name not in align.SubstitutionMatrix.list_db():
        warnings.warn("Substitution Matrix could not be found. BLOSUM62 will be used.")
        matrix = align.SubstitutionMatrix.std_protein_matrix()
    else:
        alph = seq.ProteinSequence.alphabet
        matrix = align.SubstitutionMatrix(alph, alph, name)

    return matrix

