import uuid

import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align
import biotite.sequence.align as align

from .base import BaseAligner
from ..sequences import needleman_wunsch, smith_waterman, matrices
from ..utils import enter_temp_directory


class MatchMakerAligner(BaseAligner):

    """
    Factory to configure an aligner based on
    UCSF Chimera's MatchMaker algorithms.

    Roughly, the algorithm follows these steps:

    1. Sequence alignment -- using "biotite" and "MDAnalysis"
    2. Trace (atom pairing) -- using "biotite" and "MDAnalysis"
    3. Structural superposition -- using "MDAnalysis"

    Parameters
    ----------
    alignment_strategy : str, optional, default=global
        What type of algorithm will be used to calculate the
        sequence alignment. Choose between:
        - "global" (Needleman-Wunsch)
        - "local" (Smith-Waterman)
    alignment_matrix : str, optional, default=BLOSUM62
        The substitution matrix used for scoring
    alignment_gap : int or (tuple, dtype=int), optional
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple is the gap opening penalty,
        the second integer is the gap extension penalty. The values need to be negative.
        (Default: -10)
    strict_superposition: bool, optional, default=False
        True: Will raise SelectionError if a single atom does not match between the two selections.
        False: Will try to prepare a matching selection by dropping residues with non-matching atoms.
    superposition_selection: str or AtomGroup or None, optional, default=None
        None: Apply to mobile.universe.atoms (i.e., all atoms in the context of the selection from mobile such as the rest of a protein, ligands and the surrounding water)
        str: Apply to mobile.select_atoms(selection-string), e.g "protein and name CA"
        AtomGroup: Apply to the arbitrary group of atoms
    superposition_weights: {“mass”, None} or array_like, optional
        choose weights. With "mass" uses masses as weights;
        None: weigh each atom equally
        If a float array of the same length as mobile is provided, use each element of the array_like as a weight for the corresponding atom in mobile.
    superposition_delta_mass_tolerance: float, optional, default=0.1
        Reject match if the atomic masses for matched atoms differ by more than tol_mass

    References
    ----------
    * <https://www.mdanalysis.org/docs/documentation_pages/analysis/align.html#MDAnalysis.analysis.align.alignto>

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

        self.alignment_matrix = alignment_matrix
        if self.alignment_matrix not in align.SubstitutionMatrix.list_db():
            raise ValueError(f"Substitution Matrix '{self.alignment_matrix }' could not be found.")
        else:
            self.alignment_matrix = matrices

        self.alignment_gap = alignment_gap
        self.strict_superposition = strict_superposition
        self.superposition_selection = superposition_selection
        self.superposition_weights = superposition_weights
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
        ref_sequence = self._retrieve_sequence(reference)
        mob_sequence = self._retrieve_sequence(mobile)
        alignment = self._align(ref_sequence, mob_sequence)

        # Retrieve trace atoms
        trace = alignment.trace
        # grab only rows where sequences are aligned
        trace = trace[~(trace == -1).any(axis=1)]

        aligned_residues_ref = ref_universe.residues[trace[:, 0]]
        aligned_residues_mob = mob_universe.residues[trace[:, 1]]

        # FIXME: Does MDA move the structure as part of the RMSD calculation?
        old_rmsd, new_rmsd = mda_align.alignto(
            aligned_residues_mob.atoms,
            aligned_residues_ref.atoms,
            strict=self.strict_superposition,
            select=self.superposition_selection,
            weights=self.superposition_weights,
            tol_mass=self.superposition_delta_mass_tolerance,
        )
        return {"rmsd": new_rmsd}

    @staticmethod
    def _retrieve_sequence(atomium_model):
        """
        Get the amino acid sequence

        Parameters
        ----------
        atomium model

        Returns
        -------
        str
            one-letter amino acid sequence
        """

        sequences = []
        for chain in atomium_model._chains.structures:
            sequences.append(chain.sequence)
        return "".join(sequences)

    @staticmethod
    def _atomium_to_mda_universe(atomium_model):
        with enter_temp_directory():
            filename = f"{uuid.uuid4()}.pdb"
            atomium_model.save(filename)
            return mda.Universe(filename)

    def _align(self, sequence_1, sequence_2):
        """
        Perform a global alignment, based on the the Needleman-Wunsch algorithm
        or a local alignment, based on the Smith-Waterman algorithm

        Parameters
        ----------
        sequence_1,sequence_2: str
            str of sequences

        Returns
        -------
        str
            an alignment of two sequences

        Examples
        --------
        >>> _align("seq1", "seq2","PAM250", (-5, -15))
        RKKSLVDIDLSSLRDP
        R-K-I-DLS-S-LRDP
        """

        return self._sequence_aligner(
            sequence_1, sequence_2, self.alignment_matrix, self.alignment_gap
        )
