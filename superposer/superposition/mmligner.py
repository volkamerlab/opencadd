"""
MMLigner (Collier et al., 2017) works by minimizing the ivalue of the alignment. The ivalue is based on
the Minimum Message Length framework (Wallace and Boulton, 1968; Wallace, 2005), a Bayesian framework for
statistical inductive inference. The ivalue represents the hypothetical minimum message length needed to transmit
the computed alignment losslessly (Shannon, 1948).
Using the ivalue measure, the algorithm creates crude-but-effective strucural alignments rapidly to act as seeds.
These seeds are iteratively refined over an Expectation-Maximization approach using the I-value criterion.
By substracting the ivalue from the null model, the statistical significance of the alignment can be computed. If the
difference is positive, the alignment is significant.

Collier, J.H., Allison, L., Lesk A.M., Stuckey, P.J., Garcia de la Banda , M., Konagurthu, A.S. (2017)
Statistical inference of protein structural alignments using information and compression.
Bioinformatics, 33(7), 1005-1013

Wallace,C.S. and Boulton,D.M. (1968) An information measure for classification.
Comput. J., 11, 185–194.

Wallace,C.S. (2005) Statistical and Inductive Inference Using MinimumMessage Length.
Information Science and Statistics. SpringerVerlag, New York, NY.

Shannon,C.E. (1948) A mathematical theory of communication.
Bell Syst.Tech. J., 27, 379–423.
"""

import sys
import subprocess
from copy import deepcopy
import logging

import numpy as np

# import atomium
import biotite.sequence.io.fasta as fasta

from .base import BaseAligner
from ..utils import enter_temp_directory

_logger = logging.getLogger(__name__)


class MMLignerAligner(BaseAligner):
    """
    Wraps MMLigner to superpose two protein structures.

    Parameters
    ----------
    executable : str
        Path to the MMLigner executable file
    """

    def __init__(self, executable=None):
        if executable is None:
            executable = "mmligner64.exe" if sys.platform.startswith("win") else "mmligner"
        self.executable = executable
        _logger.warning(
            "Current MMLigner wrappers produces accurate RMSD values but slightly shifted structures."
        )

    def _calculate(self, structures, *args, **kwargs):
        """
        Calculates the superposition of two protein structures.

        It is called by BaseAligner.calculate().

        Parameters
        ----------
        structures: [array like, array like]
            Sequences of two protein structures of same length

        Returns
        -------
        dict
            As returned by ``._parse_metadata(output)``.
            - ``superposed`` ([atomium.model, atomium.model]): the superposed models
            - ``scores`` (dict):
                - ``rmsd`` (float): RMSD value of the alignment
                - ``score`` (float): ivalue of the alignment
                - ``coverage`` (float): coverage of the alignment
            - ``metadata`` (dict):
                - ``alignment`` (biotite.alignment): computed alignment
                - ``rotation`` (array-like): 3x3 rotation matrix
                - ``translation`` (np.array): array containing the translation
                - ``quarternion`` (array-like): 4x4 quarternion matrix
        """

        with enter_temp_directory() as (cwd, tmpdir):
            sys.setrecursionlimit(
                100000
            )  # Needed because of the need of a copy of the structures.

            path1, path2 = self._edit_pdb(structures)
            output = subprocess.check_output(
                [self.executable, path1, path2, "-o", "temp", "--superpose"]
            )
            # We need access to the temporary files at parse time!
            result = self._parse_metadata(output.decode())
            copies = [deepcopy(structure) for structure in structures]
            superposed_models = self._calculate_transformed(copies, result["metadata"])
            result.update({"superposed": superposed_models})
        return result

    def _parse_metadata(self, output):
        """
        Retrieves RMSD, score and metadata from the output of the MMLigner subprocess.

        Parameters
        ----------
        output: str
            string of containing the stdout of the mmligener call

        Returns
        -------
        dict
            As returned by ``._parse_metadata(output)``.
            - ``scores`` (dict):
                - ``rmsd`` (float): RMSD value of the alignment
                - ``score`` (float): ivalue of the alignment
                - ``coverage`` (float): coverage of the alignment
            - ``metadata`` (dict):
                - ``alignment``: (biotite.alignment): computed alignment
                - ``rotation``: (array-like): 3x3 rotation matrix
                - ``translation``: (np.array): array containing the translation
                - ``quarternion``: (array-like): 4x4 quarternion matrix
        """
        lines = iter(output.splitlines())
        for line in lines:
            if line.startswith("RMSD"):
                rmsd = float(line.split()[2])
            elif line.startswith("Coverage"):
                coverage = float(line.split()[2])
            elif line.startswith("I(A & <S,T>)"):
                ivalue = float(line.split()[4])
            elif "Print Centers of Mass of moving set:" in line:
                moving_com = np.array([float(x) for x in next(lines).split()])
            elif "Print Centers of Mass of fixed set:" in line:
                fixed_com = np.array([float(x) for x in next(lines).split()])
            elif "Print Rotation matrix" in line:
                rotation = [[float(x) for x in next(lines).split()] for _ in range(3)]
            elif "Print Quaternion matrix" in line:
                quarternion = [[float(x) for x in next(lines).split()] for _ in range(4)]

        translation = fixed_com - moving_com

        alignment = fasta.FastaFile()
        alignment.read("temp__1.afasta")

        return {
            "scores": {"rmsd": rmsd, "score": ivalue, "coverage": coverage},
            "metadata": {
                "alignment": alignment,
                "rotation": rotation,
                "translation": translation,
                "quarternion": quarternion,
            },
        }

    def _parse_scoring(self, output):
        """
        Retrieves RMSD, score and ivalue from the output of the MMLigner subprocess.

        Parameters
        ----------
        output: str
            string containing the stdout of the mmligener call

        Returns
        -------
        dict
            As returned by ``._parse_scoring(output)``.
            - ``scores`` (dict):
                - ``rmsd`` (float): RMSD value of the alignment
                - ``score`` (float): ivalue of the alignment
                - ``coverage`` (float): coverage of the alignment
        """
        lines = iter(output.splitlines())
        for line in lines:
            if line.startswith("RMSD"):
                rmsd = float(line.split()[2])
            elif line.startswith("Coverage"):
                coverage = float(line.split()[2])
            elif line.startswith("I(A & <S,T>)"):
                ivalue = float(line.split()[4])

        return {
            "scores": {"rmsd": rmsd, "score": ivalue, "coverage": coverage},
        }

    def _calculate_transformed(self, structures, metadata):
        """
        Parse back output PDBs and construct updated atomium models.

        Parameters
        ----------
        structures: list of atomium.Model
            Original input structures

        Returns
        -------
        list of atomium.Model
            Input structures with updated coordinates
        """
        ref, original_mobile, *_ = structures
        translation = metadata["translation"]
        rotation = metadata["rotation"]

        atomium_translation = original_mobile.center_of_mass - ref.center_of_mass
        original_mobile.translate(*translation)
        original_mobile.transform(rotation)
        original_mobile.translate(ref.center_of_mass - original_mobile.center_of_mass)

        return ref, original_mobile

    def ivalue(self, structures, alignment):
        """
        Parse back output PDBs and construct updated atomium models.

        Parameters
        ----------
        structures: [array like, array like]
            sequences of two protein structures of same length
        alignment: biotite.alignment
            alignment of the given two sequences

        Returns
        -------
        dict
            As returned by ``._parse_scoring(output)``.
            - ``scores`` (dict):
                - ``rmsd`` (float): RMSD value of the alignment
                - ``score`` (float): ivalue of the alignment
                - ``coverage`` (float): coverage of the alignment
        """

        with enter_temp_directory() as (cwd, tmpdir):
            path1, path2 = self._edit_pdb(structures)

            fasta_file = fasta.FastaFile()

            for header, string in alignment.items():
                fasta_file[header] = string

            fasta_file.write("temp_alignment.afasta")

            self._edit_fasta("temp_alignment.afasta")

            output = subprocess.check_output(
                [self.executable, path1, path2, "--ivalue", "temp_alignment.afasta"]
            )
            # We need access to the temporary files at parse time!
            result = self._parse_scoring(output.decode())

        return result

    def _edit_pdb(self, structures, path=("structure1.pdb", "structure2.pdb")):
        """
        Method to write atomium protein models to PDBs readable by MMLigner.

        Parameters
        ----------
        structures: [array like, array like]
            two protein structures

        path: [str, str], Optional=["structure1.pdb, "structure2.pdb"]
            Path where the pdbs should be written

        Returns
        -------
        str, str
            Paths of both structures

        .. note::

            This is a temporary workaround to fix issue #9 at:
            https://github.com/volkamerlab/superposer/issues/9
        """
        assert len(path) == 2

        structures[0].save(path[0])
        structures[1].save(path[1])

        for i in range(len(path)):
            pdb = []
            with open(path[i], "r") as infile:
                pdb = infile.readlines()
                for j in range(1, len(pdb)):
                    if pdb[j].startswith("TER"):
                        pdb[j] = pdb[j].split()[0] + "    " + pdb[j - 1].split()[1] + "\n"

            self._write_pdb(path[i], pdb)

        return path[0], path[1]

    def _write_pdb(self, path, pdb):
        """
        Method to write atomium protein models to PDBs readable by MMLigner.

        Parameters
        ----------
        path: str
            Path where the pdb should be written

        pdb: array-like
            edited pdb file

        .. note::

            This is a temporary workaround to fix issue #9 at:
            https://github.com/volkamerlab/superposer/issues/9
        """
        with open(path, "w") as outfile:
            for line in pdb:
                outfile.write(line)

    def _edit_fasta(self, path):
        """
        Method to edit FASTA files written by biotite to FASTA files readable by MMLigner. This is needed,
        because MMLigner expects an empty line after each sequence.

        Parameters
        ----------
        path: str
            Path to the fasta file that is to be edited.
        """
        with open(path, "r") as fasta:
            data = fasta.readlines()

        lines = iter(data)
        for line in lines:
            if next(lines).startswith(">structure2.pdb"):
                line = line + "\n"

        data[-1] = data[-1] + "\n"

        with open(path, "w") as fasta:
            for line in data:
                fasta.write(line)
