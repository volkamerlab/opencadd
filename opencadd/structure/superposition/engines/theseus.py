"""
Theseus superposes a set of macromolecular structures simultaneously
using the method of maximum like-lihood (ML), rather than the
conventional least-squares criterion. Theseus assumes that the
structures are distributed according to a matrix Gaussian
distribution and that the eigenvalues of the atomic
covariancematrix are hierarchically distributed according
to an inverse gamma distribution. This ML superpositioning model
produces much more accurate results by essentially downweighting
variable regions of the structuresand by correcting for correlations
among atoms.

References
----------
* https://theobald.brandeis.edu/theseus/
* Optimal simultaneous superpositioning of multiple structures with missing data.
  Theobald, Douglas L. & Steindel, Philip A. (2012) Bioinformatics 28 (15): 1972-1979
* Accurate structural correlations from maximum likelihood superpositions.
  Theobald, Douglas L. & Wuttke, Deborah S. (2008) PLOS Computational Biology 4(2):e43
"""
import os
import subprocess
import logging
from distutils.spawn import find_executable

import numpy as np

from .base import BaseAligner
from ....utils import enter_temp_directory


_logger = logging.getLogger(__name__)


class TheseusAligner(BaseAligner):

    """
    Superpose structures with different sequences but similar structures

    Parameters
    ----------
    alignment_max_iterations : int
        number of iterations for alignment program (muscle)
    statistics_only : bool
        if True, add `-I` flag to just compute the statistics (no superposition)
    """

    def __init__(self, alignment_max_iterations: int = 32, statistics_only: bool = False) -> None:
        self.alignment_max_iterations = alignment_max_iterations
        self.statistics_only = statistics_only

        self._fastafile = "theseus.fasta"
        self._filemap_file = "theseus.filemap"
        self._alignment_file = "theseus.aln"
        self._alignment_log = "muscle.log"
        self._alignment_file_biotite = "theseus_biotite.aln"
        self._alignment_executable = "muscle"
        self._theseus_transformation_file = "theseus_transf.txt"

    def _safety_checks(self):
        """
        Check if `theseus` is installed (executable found?).

        Raises
        ------
        OSError
            Raises error if executable `theseus` cannot be found.
        """

        theseus = find_executable("theseus")
        if theseus is None:
            raise OSError("theseus cannot be located. Is it installed?")
        # proceed normally

    def _calculate(self, structures, *args, **kwargs) -> dict:
        """
        Align the sequences with an alignment tool (``muscle``)

        A mode for superpositioning macromolecules with
        identical sequences and numbers of residues,
        for instance, multiple models in an NMR family
        or multiple structures from different crystal
        forms of the same protein.


        Parameters
        ----------
        structures : list
            list of opencadd.core.Structures objects
        **kwargs : dict
            optional parameters

        Returns
        -------
        dict
            As returned by ``._parse_superposition(output)``.

            - ``scores`` (dict)
                    ``rmsd`` (float): RMSD value of the alignment
            - ``metadata`` (dict)
                    ``least_squares``: least squares
                    ``maximum_likelihood``: maximum likelihood
                    ``log_marginal_likelihood``: log marginal likelihood
                    ``aic``: Akaike information criterion
                    ``bic``: Bayesion information criterion
                    ``omnibus_chi_square``: omnibus chi square
                    ``hierarchical_var_chi_square``: hierarchical var chi square
                    ``rotational_translational_covar_chi_square``: rotational translational covar chi square
                    ``hierarchical_minimum_var``: hierarchical minimum var
                    ``hierarchical_minimum_sigma``: hierarchical minimum sigma
                    ``skewness``: skewness
                    ``skewness_z``: skewness z
                    ``kurtosis``: kurtosis
                    ``kurtosis_z``: kurtosis z
                    ``data_pts``: data points
                    ``free_params``: free params
                    ``d_p``: d_p
                    ``median_structure``: median_structure
                    ``n_total``: number total
                    ``n_atoms``: number of atoms
                    ``n_structures``: number of structures
                    ``total_rounds``: total rounds
        """

        with enter_temp_directory(remove=True) as (cwd, tmpdir):
            _logger.debug("All files are located in: %s", tmpdir)

            # 1st - Dump PDBs to disk
            filenames = []
            for index, structure in enumerate(structures):
                filename = f"structure_{index}.pdb"
                structure.write(filename)
                filenames.append(filename)

            # 2nd - Prepare sequence alignment with MUSCLE
            seq_alignment_output = self._run_alignment(filenames)
            _logger.info(seq_alignment_output)

            # 3rd - Run theseus superposition itself
            theseus_output = subprocess.check_output(
                [
                    "theseus",
                    "-f",
                ]
                + (["-I"] if self.statistics_only else [])
                + [
                    "-M",
                    self._filemap_file,
                    "-A",
                    self._alignment_file,
                    "-o",  # all matrices with respect to this file
                    filenames[0],
                    *filenames,
                ],
                universal_newlines=True,
            )
            _logger.debug(theseus_output)
            results = self._parse_superposition(theseus_output)

        # 4th reapply transformation to the mobile model
        mobile = structures[1]
        transformation = results["metadata"]["transformation"]
        rotation = transformation[:, :3]
        translation = transformation[:, 3:].reshape((3,))
        mobile.atoms.translate(translation)
        mobile.atoms.rotate(rotation)
        results["superposed"] = structures
        return results

    def _get_fasta(self, filenames) -> str:
        """
        Use Theseus to create fasta files

        Parameter
        ---------
        filenames : list
            list of pdb filenames

        Returns
        -------
        str
            Output of the created fasta filenames
        """
        return subprocess.check_output(["theseus", "-f", "-F", *filenames])

    def _concatenate_fasta(self, filenames) -> None:
        """
        Concatenate the fasta files created with Theseus into one multiple sequence fasta file

        Parameter
        ---------
        filenames : list
            list of pdb filenames
        """
        with open(self._fastafile, "w") as outfile:
            for fname in filenames:
                with open(f"{fname}.fst") as infile:
                    for line in infile:
                        outfile.write(line)

    def _filemap(self, filenames) -> None:
        """
        Create filemap for Theseus

        Parameter
        ---------
        filenames : list
            list of pdb filenames
        """
        with open(self._filemap_file, "w") as outfile:
            for fname in filenames:
                outfile.write(f"{fname} {fname}\n")

    def _run_alignment(self, filenames):
        """
        Run MUSCLE

        Returns
        -------
        filenames : list of str
            Paths to PDB files containing the structures
        str
            Output of MUSCLE
        """
        self._get_fasta(filenames)
        self._concatenate_fasta(filenames)
        self._filemap(filenames)
        output = subprocess.check_output(
            [
                self._alignment_executable,
                "-maxiters",
                str(self.alignment_max_iterations),
                "-in",
                self._fastafile,
                "-out",
                self._alignment_file,
                "-clwstrict",
                "-verbose",
                "-log",
                self._alignment_log,
            ],
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        _logger.info(output)
        return output

    def _get_transformation_matrix(self) -> dict:
        """
        Get the rotation matrix and translation vector

        Returns
        -------
        dict
            Rotation matrix and translation vector
        """
        translations, rotations = {}, {}
        if not os.path.isfile(self._theseus_transformation_file):
            return None
        with open(self._theseus_transformation_file, "r") as infile:
            lines = infile.readlines()
            for line in lines:
                if " R:" in line:  # rotation matrix
                    model, matrix = line.split(":")
                    model_id = int(model.split()[1])
                    matrix_array = np.reshape(list(map(float, matrix.split())), (-1, 3))
                    rotations[model_id] = matrix_array
                elif " t:" in line:  # translation vector
                    model, vector = line.split(":")
                    model_id = int(model.split()[1])
                    vector_array = np.array([[float(x)] for x in vector.split()])
                    translations[model_id] = vector_array
        matrices = {}
        for model_id, rotation in rotations.items():
            translation = translations[model_id]
            matrix = np.empty((3, 4))
            matrix[:, :3] = rotation
            matrix[:, 3:] = translation
            matrices[model_id] = matrix
        return matrices[2]

    def _parse_superposition(self, output) -> dict:
        """
        Parse the output from theseus itself

        Parameters
        ----------
        output : bytes

        Returns
        -------
        dict
            All the information provided by Theseus
        """
        for line in output.splitlines():
            blocks = line.split()
            if "Classical LS pairwise <RMSD>" in line:
                rmsd = float(blocks[5])
            elif "Least-squares <sigma>" in line:
                least_squares = float(blocks[3])
            elif "Maximum Likelihood <sigma>" in line:
                maximum_likelihood = float(blocks[4])
            elif "Marginal Log Likelihood" in line:
                log_marginal_likelihood = float(blocks[4])
            elif "AIC" in line:
                aic = float(blocks[2])
            elif "BIC" in line:
                bic = float(blocks[2])
            elif "Omnibus chi^2" in line:
                omnibus_chi_square = float(blocks[3])
            elif "Hierarchical var" in line:
                hierarchical_var_chi_square = float(blocks[6])
            elif "Rotational, translational, covar" in line:
                rotational_translational_covar_chi_square = float(blocks[5])
            elif "Hierarchical minimum var" in line:
                hierarchical_minimum_var = float(blocks[5])
                hierarchical_minimum_sigma = float(blocks[-1][1:-1])
            elif "skewness Z-value" in line:
                skewness_z = float(blocks[3])
            elif "skewness" in line:
                skewness = float(blocks[2])
            elif "kurtosis Z-value" in line:
                kurtosis_z = float(blocks[3])
            elif "kurtosis" in line:
                kurtosis = float(blocks[2])
            elif "Data pts" in line:
                fields = line.split(",")
                data_pts = float(fields[0].split()[-1])
                free_params = float(fields[1].split()[-1])
                d_p = float(fields[2].split()[-1])
            elif "Median" in line:
                median_structure = float(blocks[4][1:])
            elif "N(total)" in line:
                fields = line.split(",")
                n_total = float(fields[0].split()[-1])
                n_atoms = float(fields[1].split()[-1])
                n_structures = float(fields[2].split()[-1])
            elif "Total rounds" in line:
                total_rounds = float(blocks[3])
        return {
            "scores": {"rmsd": rmsd},
            "metadata": {
                "transformation": self._get_transformation_matrix(),
                "least_squares": least_squares,
                "maximum_likelihood": maximum_likelihood,
                "log_marginal_likelihood": log_marginal_likelihood,
                "aic": aic,
                "bic": bic,
                "omnibus_chi_square": omnibus_chi_square,
                "hierarchical_var_chi_square": hierarchical_var_chi_square,
                "rotational_translational_covar_chi_square": rotational_translational_covar_chi_square,
                "hierarchical_minimum_var": hierarchical_minimum_var,
                "hierarchical_minimum_sigma": hierarchical_minimum_sigma,
                "skewness": skewness,
                "skewness_z": skewness_z,
                "kurtosis": kurtosis,
                "kurtosis_z": kurtosis_z,
                "data_pts": data_pts,
                "free_params": free_params,
                "d_p": d_p,
                "median_structure": median_structure,
                "n_total": n_total,
                "n_atoms": n_atoms,
                "n_structures": n_structures,
                "total_rounds": total_rounds,
            },
        }
