"""
opencadd.databases.klifs.exceptions

Defines KLIFS exceptions.
"""


class KlifsPocketIncompleteError(Exception):
    """
    Exception raised for errors in the KLIFS pocket length.

    Attributes
    ----------
    pocket_length : int
        Pocket length that shall raise an error.
    message : str
        Explanation of the error.
    """

    def __init__(self, pocket_length, message="Length of KLIFS pocket sequence but must be 85"):
        self.pocket_length = pocket_length
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} but is {self.pocket_length}"


class KlifsPocketUnequalSequenceStructure(Exception):
    """
    Exception raised for errors in case the KLIFS pocket length is unequal in the sequence and the
    structure (file).

    Attributes
    ----------
    pocket_length : int
        Pocket length that shall raise an error.
    message : str
        Explanation of the error.
    """

    def __init__(
        self,
        pocket_length_sequence,
        pocket_length_structure,
        message="Length of KLIFS pocket is unequal in sequence and structure",
    ):
        self.pocket_length_sequence = pocket_length_sequence
        self.pocket_length_structure = pocket_length_structure
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return (
            f"{self.message}: "
            f"Sequence has {self.pocket_length_sequence} and "
            f"structure has {self.pocket_length_structure} residues."
        )
