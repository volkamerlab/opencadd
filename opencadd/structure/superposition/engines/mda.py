"""
Aligner based on MDAnalysis' superposition algorithms.
"""

import logging

import numpy as np
from MDAnalysis.analysis import align as mda_align, rms
from MDAnalysis.lib.util import canonical_inverse_aa_codes, convert_aa_code
import biotite.sequence.align as align
from biotite.sequence.io.fasta import FastaFile

from .base import BaseAligner
from ..sequences import sequence_alignment, fasta2select
from ....utils import enter_temp_directory

_logger = logging.getLogger(__name__)


class MDAnalysisAligner(BaseAligner):

    """
    Factory to configure an aligner based on
    MDAnalysis' superposition algorithms.

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
    alignment_gap : int or (tuple, dtype=int), optional, default=-10
        Int the value will be interpreted as general gap penalty.
        Tupel is provided, an affine gap penalty is used. The first integer in the tuple
        is the gap opening penalty, the second integer is the gap extension penalty.
        The values need to be negative.
    strict_superposition: bool, optional, default=False
        True: Will raise SelectionError if a single atom does not match between the two selections.
        False: Will try to prepare a matching selection by dropping residues with non-matching atoms.
    per_residue_selection: str or AtomGroup or None, optional, default=None
        None: Apply to mobile.universe.atoms (i.e., all atoms in the context of the selection from
        mobile such as the rest of a protein, ligands and the surrounding water)
        str: Apply to mobile.select_atoms(selection-string), e.g "protein and name CA"
        AtomGroup: Apply to the arbitrary group of atoms
    superposition_weights: {“mass”, None} or array_like, optional
        choose weights. With "mass" uses masses as weights;
        None: weigh each atom equally
        If a float array of the same length as mobile is provided, use each element of the
        array_like as a weight for the corresponding atom in mobile.
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
        per_residue_selection="name CA and not altloc B and not altloc C",  # TODO: why not "name CA and altloc A"
        superposition_weights=None,
        superposition_delta_mass_tolerance=0.1,
    ):

        self.alignment_strategy = alignment_strategy.lower()
        if self.alignment_strategy == "global":
            self._align_local = False
        elif self.alignment_strategy == "local":
            self._align_local = True
        else:
            raise ValueError("`alignment_strategy` must be one of `global, local`.")

        if alignment_matrix not in align.SubstitutionMatrix.list_db():
            raise ValueError(f"Substitution Matrix '{alignment_matrix }' could not be found.")
        else:
            self.alignment_matrix = alignment_matrix

        self.alignment_gap = alignment_gap
        self.strict_superposition = strict_superposition
        self.per_residue_selection = per_residue_selection
        self.superposition_weights = superposition_weights
        self.superposition_delta_mass_tolerance = superposition_delta_mass_tolerance

    def _safety_checks(self):
        """
        Check for `mda` installation passes; added here only for consistency across the engines.
        """

        pass

    def _calculate(self, structures, selections, *args, **kwargs):
        """

        Parameters
        ----------
        structures : list of opencadd.core.Structure
            First one will be the target (static structure). Following ones will be mobile.
        selections : list of MDAnalysis.core.groups.AtomGroup
            Selection of atoms on which the calculation is done. If no selection is given, the whole structure is used for calculation.

        Returns
        -------
        dict
            superposed models
            rmsd
            metadata
        """

        ref_universe, mob_universe = structures

        # we want to perform the calculation on the selections
        # if there is no selection, we reference to the structures
        # e.g. id(selections[0]) equals id(structures[0])
        # after calculation is done, we perform the transformation on the structure
        # this way, we do not throw away any atoms and transform the whole structure
        if not selections:
            selections.append(ref_universe)
            selections.append(mob_universe)

        ref_selection, mob_selection = selections

        # Get matching atoms of selection
        selection, alignment = self.matching_selection(*selections)
        ref_atoms = ref_selection.select_atoms(selection["reference"])
        mobile_atoms = mob_selection.select_atoms(selection["mobile"])
        if len(ref_atoms) == len(mobile_atoms):
            coverage = len(ref_atoms)
        else:
            raise ValueError(
                "The number of atoms to match should be the same for both structures."
            )

        # Compute initial RMSD (no preprocessing)
        initial_rmsd = rms.rmsd(ref_atoms.positions, mobile_atoms.positions)

        # Compute centered RMSD (both structures share now the same center of mass)
        ref_weights = np.asarray([a.mass for a in ref_atoms])
        mobile_weights = np.asarray([a.mass for a in mobile_atoms])
        ref_com = ref_atoms.center(ref_weights)
        mobile_com = mobile_atoms.center(mobile_weights)
        ref_coordinates = ref_atoms.positions - ref_com
        mobile_coordinates = mobile_atoms.positions - mobile_com
        centered_rmsd = rms.rmsd(ref_coordinates, mobile_coordinates)

        # Calculate optimum rotation matrix
        rotation, rmsd = mda_align.rotation_matrix(mobile_coordinates, ref_coordinates)
        # Transformation on the whole structure
        mob_universe.atoms.translate(-mobile_com)
        mob_universe.atoms.rotate(rotation)

        return {
            "superposed": [ref_universe, mob_universe],
            "scores": {"rmsd": rmsd, "coverage": coverage},
            "metadata": {
                "selection": selection,
                "alignment": alignment,
                "initial_rmsd": initial_rmsd,
                "centered_rmsd": centered_rmsd,
                "translation": ref_com,
                "rotation": rotation,
            },
        }

    def matching_selection(self, ref_selection, mob_selection):
        """
        Compute best matching atom sets

        Parameters
        ----------
        selections : list of MDAnalysis.core.groups.AtomGroup

        Returns
        -------
        dict
            Two-element dictionary with the selection string
            to obtain the matching atoms on the original structures
        alignment : biotite.Alignment
            The sequence alignment
        """
        # Compute sequence alignment and matching atoms
        ref_sequence, ref_resids, ref_segids = self._retrieve_sequence(ref_selection)
        mob_sequence, mob_resids, mob_segids = self._retrieve_sequence(mob_selection)
        alignment = self._align(ref_sequence, mob_sequence)
        with enter_temp_directory():
            fasta = FastaFile()
            fasta["ref"], fasta["mob"], *_empty = alignment.get_gapped_sequences()
            fasta.write("temp.fasta")
            selection = fasta2select(
                "temp.fasta",
                ref_resids=ref_resids,
                target_resids=mob_resids,
                ref_segids=ref_segids,
                target_segids=mob_segids,
                backbone_selection=self.per_residue_selection,
            )
        return selection, alignment

    @staticmethod
    def _retrieve_sequence(atom_group):
        """
        Get the amino acid sequence

        Parameters
        ----------
        atom_group : MDAnalysis.core.groups.AtomGroup

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
        backbone_atom_names = {"C", "CA", "N", "O"}
        incomplete_residues = []
        for segment in atom_group.segments:
            for residue in segment.residues:
                if residue.resname in canonical_inverse_aa_codes:
                    residue_atom_names = set([a.name for a in residue.atoms])
                    if not backbone_atom_names.issubset(residue_atom_names):
                        incomplete_residues.append(residue)
                    # some structures are missing atoms in the aligned residues,
                    # then an error would occur, where the amount of atoms of both
                    #  structures is not the same (needs to be the same amount for superposition)
                    else:
                        sequences.append(convert_aa_code(residue.resname))
                        residue_ids.append(residue.resid)
                        segment_ids.append(residue.segid)
        if incomplete_residues:
            _logger.warning(
                "%d residues in %s are missing backbone atoms. "
                "If this system was obtained from a larger structure using a "
                "selection, consider wrapping such selection with "
                "`same residue as (<your original selection>)` "
                "to avoid potential matching problems.",
                len(incomplete_residues),
                atom_group,
            )
        return "".join(sequences), residue_ids, segment_ids

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

        return sequence_alignment(
            sequence_1,
            sequence_2,
            self.alignment_matrix,
            self.alignment_gap,
            local=self._align_local,
        )
