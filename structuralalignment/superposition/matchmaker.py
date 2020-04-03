import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align

from .base import BaseAligner
from ..sequences import needleman_wunsch, smith_waterman


class MatchMakerAligner(BaseAligner):
    """
    Factory to configure an aligner based on
    UCSF Chimera's MatchMaker algorithms.

    Roughly, the algorithm follows these steps:

    1. Sequence alignment -- using ``biotite``
    2. Trace (atom pairing) -- using ``biotite`` and ``MDAnalysis``
    3. Structural superposition -- using ``MDAnalysis``

    Parameters
    ----------
    alignment_strategy : str, optional, default=global
        What type of algorithm will be used to calculate the
        sequence alignment. Choose between:
        - ``global`` (Needleman-Wunsch)
        - ``local`` (Smith-Waterman)
    alignment_matrix : str, optional, default=BLOSUM62
        The substitution matrix used for scoring
    alignment_gap : int or (tuple, dtype=int), optional
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.
        (Default: -10)
    strict_superposition
        TODO: add description and typing
    superposition_selection
        TODO: add description and typing
    superposition_weights
        TODO: add description and typing
    superposition_delta_mass_tolerance
        TODO: add description and typing

    References
    ----------
    * <LINK>!
    """

    def __init__(
        self,
        alignment_strategy: str = "global",
        alignment_matrix: str = "BLOSUM62",
        alignment_gap: int = -10,
        strict_superposition: bool = False,
        superposition_selection="protein and name CA",
        superposition_weights=None,
        superposition_delta_mass_tolerance=0.1,
    ):
        self.alignment_strategy = alignment_strategy.lower()
        if self.alignment_strategy == "global":
            self._sequence_aligner = needleman_wunsch
        elif self.alignment_strategy == "local":
            self._sequence_aligner = smith_waterman
        else:
            raise ValueError("`alignment_strategy` must be one of `global, local`.")

        # TODO: Add checks for allowed values (as above) and raise ValueError with informative message otherwise
        self.alignment_matrix = alignment_matrix
        # TODO: Add checks for allowed values (as above) and raise ValueError with informative message otherwise
        self.alignment_gap = alignment_gap
        # TODO: Add checks for allowed values (as above) and raise ValueError with informative message otherwise
        self.strict_superposition = strict_superposition
        # TODO: Add checks for allowed values (as above) and raise ValueError with informative message otherwise
        self.superposition_selection = superposition_selection
        # TODO: Add checks for allowed values (as above) and raise ValueError with informative message otherwise
        self.superposition_weights = superposition_weights
        # TODO: Add checks for allowed values (as above) and raise ValueError with informative message otherwise
        self.superposition_delta_mass_tolerance = superposition_delta_mass_tolerance

    def _calculate(self, structures):
        """

        Parameters
        ----------
        structures : list of atomium.Model
            First one will be the target (static structure). Following, will be mobile.

        Returns
        -------
        dict
            superposed models
            rmsd
            metadata
        """
        if len(structures) > 2:
            raise NotImplementedError(
                "This method can only be used for two structures at the same time, for now"
            )
        reference, mobile = structures
        ref_universe = self._atomium_to_mda_universe(reference)
        mob_universe = self._atomium_to_mda_universe(mobile)

        # Compute sequence alignment
        ref_sequence = _retrieve_sequence(reference)
        mob_sequence = _retrieve_sequence(mobile)
        alignment = self._align(ref_sequence, mob_sequence)

        # Retrieve trace atoms
        trace = alignment.trace
        # Filter residue pairs that are == -1, which means they are a gap (not aligned!)
        trace = trace[~(trace == -1).any(axis=1)]

        aligned_residues_ref = ref_universe.residues[trace[:, 0]]
        aligned_residues_mob = mob_universe.residues[trace[:, 1]]

        # FIXME: Does MDA move the structure as part of the RMSD calculation?
        old_rmsd, new_rmsd = mda_align.alignto(
            aligned_residues_mob.atoms,
            aligned_residues_ref.atoms,
            strict=self.strict_superposition,
            select=self.superposition_selection,
        )

    @staticmethod
    def _retrieve_sequence(atomium_model):
        sequences = []
        for chain in atomium_model._chains.structures:
            sequences.append(chain.sequence)
        return "".join(sequences)

    @staticmethod
    def _atomium_to_mda_universe(atomium_model):
        # FIXME: This is not implemented yet!
        raise NotImplementedError("TODO: Jaime will take care of this")
        return mda.Universe(...)

    def _align(self, sequence_1, sequence_2):
        """
        Examples
        --------
        >>> get_alignment("xxxx.pdb", "xxxx.pdb", False, "PAM250", (-5, -15))

        RKKSLVDIDLSSLRDP
        R-K-I-DLS-S-LRDP
        """
        return self._sequence_aligner(
            sequence_1, sequence_2, self.alignment_matrix, self.alignment_gap
        )
