"""
Aligner based on UCSF Chimera's MatchMaker algorithms.
"""

import uuid
from copy import deepcopy

import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align, rms
from MDAnalysis.lib.util import canonical_inverse_aa_codes, convert_aa_code
import biotite.sequence.align as align
from biotite.sequence.io.fasta import FastaFile

from .base import BaseAligner
from ..sequences import needleman_wunsch, smith_waterman, fasta2select
from ..utils import enter_temp_directory


class MatchMakerAligner(BaseAligner):

    """
    Factory to configure an aligner based on
    UCSF Chimera's MatchMaker algorithms.

    Roughly, the algorithm follows these steps:

    1. Sequence alignment -- using "biotite" and "MDAnalysis"
    2. Structural superposition -- using "MDAnalysis"

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
        superposition_selection="name CA",
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

        if alignment_matrix not in align.SubstitutionMatrix.list_db():
            raise ValueError(
                f"Substitution Matrix '{alignment_matrix }' could not be found."
            )
        else:
            self.alignment_matrix = alignment_matrix

        self.alignment_gap = alignment_gap
        self.strict_superposition = strict_superposition
        self.superposition_selection = superposition_selection
        self.superposition_weights = superposition_weights
        self.superposition_delta_mass_tolerance = superposition_delta_mass_tolerance

    # pylint: disable=arguments-differ
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
        mob_universe_cp = self._atomium_to_mda_universe(mobile)

        # Compute sequence alignment
        ref_sequence, ref_resids, ref_segids = self._retrieve_sequence(ref_universe)
        mob_sequence, mob_resids, mob_segids = self._retrieve_sequence(mob_universe)
        alignment = self._align(ref_sequence, mob_sequence)
        with enter_temp_directory():
            fasta = FastaFile()
            fasta["ref"] = alignment.get_gapped_sequences()[0]
            fasta["mob"] = alignment.get_gapped_sequences()[1]
            fasta.write("temp.fasta")
            selection = fasta2select(
                "temp.fasta",
                is_aligned=True,
                ref_resids=ref_resids,
                target_resids=mob_resids,
                ref_segids=ref_segids,
                target_segids=mob_segids,
                backbone_selection=self.superposition_selection,
            )

        ref_atoms = ref_universe.select_atoms(selection["reference"])
        mobile_atoms = mob_universe.select_atoms(selection["mobile"])
        initial_rmsd = rms.rmsd(ref_atoms.positions, mobile_atoms.positions)

        mobile_atoms.translate(-mobile_atoms.center_of_mass())
        ref_atoms.translate(-ref_atoms.center_of_mass())
        rotation, rmsd = mda_align.rotation_matrix(
            ref_atoms.positions, mobile_atoms.positions
        )
        mob_universe.atoms.translate(-mob_universe.atoms.center_of_mass())
        mob_universe.atoms.rotate(rotation)
        mob_universe.atoms.translate(ref_universe.atoms.center_of_mass())

        return {
            "superposed": [ref_universe, mob_universe, mob_universe_cp],
            "scores": {"rmsd": rmsd},
            "metadata": {
                "selection": selection,
                "alignment": alignment,
                "initial_rmsd": initial_rmsd,
            },
        }

    @staticmethod
    def _retrieve_sequence(universe):
        """
        Get the amino acid sequence

        Parameters
        ----------
        universe : mdanalysis.Universe

        Returns
        -------
        str
            one-letter amino acid sequence
        tuple of int
            residue ids of the protein sequence
        """
        sequences = []
        residue_ids = []
        segment_ids = []
        protein = universe.select_atoms("protein")
        for segment in protein.segments:
            for residue in segment.residues:
                if residue.resname in canonical_inverse_aa_codes:
                    sequences.append(convert_aa_code(residue.resname))
                    residue_ids.append(residue.resid)
                    segment_ids.append(residue.segid)
        return "".join(sequences), residue_ids, segment_ids

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
