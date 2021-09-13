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
import logging
from distutils.spawn import find_executable

import numpy as np
import biotite.sequence.io.fasta as fasta

from .base import BaseAligner
from ....utils import enter_temp_directory

_logger = logging.getLogger(__name__)


class MMLignerAligner(BaseAligner):
    """
    Wraps MMLigner to superpose two protein structures.

    Parameters
    ----------
    executable : str
        Path to the MMLigner executable file
    protein_selection : str, optional, default=protein
        MMLigner will not accept residues beyond the 20 standard AA.

    """

    def __init__(self, executable=None, protein_selector="protein"):
        if executable is None:
            executable = "mmligner64.exe" if sys.platform.startswith("win") else "mmligner"
        self.executable = executable
        self.protein_selector = protein_selector
        _logger.warning(
            "Current MMLigner wrappers produces accurate RMSD values but slightly shifted structures!"
        )

    def _safety_checks(self):
        """
        Check if `mmligner` is installed (executable found?).

        Raises
        ------
        OSError
            Raises error if executable `mmligner` cannot be found.
        """

        mmligner = find_executable("mmligner")
        if mmligner is None:
            raise OSError("mmligner cannot be located. Is it installed?")
        # proceed normally

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
            - ``superposed`` ([opencadd.core.Structure, opencadd.core.Structure]): the superposed models
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
            # Needed because of the need of a copy of the structures.
            sys.setrecursionlimit(100000)

            path1, path2 = self._edit_pdb(structures)
            output = subprocess.check_output(
                [self.executable, path1, path2, "-o", "temp", "--superpose"]
            )
            # We need access to the temporary files at parse time!
            result = self._parse_metadata(output.decode())

            # checks if there is metadata in the dict, if not, there was no significant alignment found.
            if "metadata" in result:
                superposed_models = self._calculate_transformed(structures, result["metadata"])
                result["superposed"] = superposed_models
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
                quaternion = [[float(x) for x in next(lines).split()] for _ in range(4)]

        # checks if there is a signifcant alignment
        if rmsd == 0 and coverage == 0:
            return {"scores": {"rmsd": rmsd, "score": ivalue, "coverage": coverage}}
        else:
            # fixed_com, moving_com, rotation and quaternion can only be obtained
            # if the patched mmligner is used (check /devtools/conda-recipes/mmligner)
            # -- this will fail in CI for now --
            translation = fixed_com - moving_com

            alignment = fasta.FastaFile()
            alignment.read("temp__1.afasta")

            return {
                "scores": {"rmsd": rmsd, "score": ivalue, "coverage": coverage},
                "metadata": {
                    "alignment": alignment,
                    "rotation": rotation,
                    "translation": translation,
                    "quaternion": quaternion,
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
        Parse back output PDBs and construct updated Structure objects.

        Parameters
        ----------
        structures: list of opencadd.core.Structure
            Original input structures

        Returns
        -------
        list of opencadd.core.Structure
            Input structures with updated coordinates
        """
        ref, mobile, *_ = structures
        translation = metadata["translation"]
        rotation = metadata["rotation"]

        mob_com = mobile.atoms.center_of_geometry()
        ref_com = ref.atoms.center_of_geometry()

        mobile.atoms.translate(-mob_com)
        mobile.atoms.rotate(rotation)
        mobile.atoms.translate(ref_com)

        return ref, mobile

    def ivalue(self, structures, alignment):
        """
        Parse back output PDBs and construct updated Structure models.

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
            paths = "structure1.pdb", "structure2.pdb"
            structures[0].select_atoms(self.protein_selector).write(paths[0])
            structures[1].select_atoms(self.protein_selector).write(paths[1])

            fasta_file = fasta.FastaFile()

            for header, string in alignment.items():
                fasta_file[header] = string

            fasta_file.write("temp_alignment.afasta")

            self._edit_fasta("temp_alignment.afasta")

            output = subprocess.check_output(
                [self.executable, paths[0], paths[1], "--ivalue", "temp_alignment.afasta"]
            )
            # We need access to the temporary files at parse time!
            result = self._parse_scoring(output.decode())

        return result

    def _edit_pdb(self, structures, path=("structure1.pdb", "structure2.pdb")):
        """
        Method to write Structure protein models to PDBs readable by MMLigner.

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
            https://github.com/volkamerlab/opencadd/issues/9
        """
        assert len(path) == 2

        structures[0].select_atoms(self.protein_selector).write(path[0])
        structures[1].select_atoms(self.protein_selector).write(path[1])

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
        Method to write Structure protein models to PDBs readable by MMLigner.

        Parameters
        ----------
        path: str
            Path where the pdb should be written

        pdb: array-like
            edited pdb file

        .. note::

            This is a temporary workaround to fix issue #9 at:
            https://github.com/volkamerlab/opencadd/issues/9
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
