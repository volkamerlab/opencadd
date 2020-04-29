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

import subprocess
import logging
import atomium

import numpy as np

from structuralalignment.superposition.base import BaseAligner
from structuralalignment.utils import enter_temp_directory


_logger = logging.getLogger(__name__)


class TheseusAligner(BaseAligner):

    """
    Superpose structures with different sequences but similar structures

    Parameters
    ----------
    max_iterations : int
        number of iterations for muscle
    """

    def __init__(self, max_iterations: int = 32) -> None:
        self.max_iterations = max_iterations

        self.fastafile = "theseus.fasta"
        self.filemap_file = "theseus.filemap"
        self.alignment_file = "theseus.aln"
        self.alignment_file_biotite = "theseus_biotite.aln"
        self.alignment_app = "muscle"

    def _create_pdb(self, structures) -> list:
        """
        Create files with atomium models

        Parameter
        ---------
        structures : list
            list of atomium models

        Returns
        -------
        list
            list of pdb filenames
        """
        fetched_pdbs_filename = []
        for index, pdbs in enumerate(structures):
            pdb_filename = f"{index}.pdb"
            pdbs.save(pdb_filename)
            fetched_pdbs_filename.append(pdb_filename)
        return fetched_pdbs_filename

    def _get_fasta(self, pdbs_filename) -> str:
        """
        Use Theseus to create fasta files

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames

        Returns
        -------
        str
            Output of the created fasta filenames
        """
        return subprocess.check_output(["theseus", "-f", "-F", *pdbs_filename])

    def _concatenate_fasta(self, pdbs_filename) -> None:
        """
        Concatenate the fasta files created with Theseus into one multiple sequence fasta file

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames
        """
        with open(self.fastafile, "w") as outfile:
            for fname in pdbs_filename:
                with open(f"{fname}.fst") as infile:
                    for line in infile:
                        outfile.write(line)

    def _filemap(self, pdbs_filename) -> None:
        """
        Create filemap for Theseus

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames
        """
        with open(self.filemap_file, "w") as outfile:
            for fname in pdbs_filename:
                outfile.write(f"{fname} {fname}\n")

    # pylint: disable=arguments-differ
    def _calculate(self, structures, identical: bool = False, **kwargs) -> dict:
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
            list of atomium objects
        identical : bool, optional, default=False
            Whether to use sequence alignment to obtain an optimum pairing for
            different length sequences or not (recommended: assume sequences are different).
        **kwargs : dict
            optional parameter

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

        with enter_temp_directory(remove=False) as (cwd, tmpdir):
            _logger.info("All files are located in: %s", tmpdir)

            pdbs_filename = self._create_pdb(structures)

            if identical:
                superposition_output = self._run_theseus_identical(pdbs_filename)
            else:
                superposition_output = self._run_theseus_different(pdbs_filename)
                _logger.info(superposition_output)
            superposed_pdb_models = self._get_superposed_models(pdbs_filename)
            transformation_matrix = self._get_transformation_matrix()
        results = self._parse_superposition(superposition_output)
        results["superposed"] = superposed_pdb_models
        results["metadata"]["transformation"] = transformation_matrix
        return results

    def _run_theseus_identical(self, pdbs_filename) -> str:
        """
        Superpose identical sequences with Theseus

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames

        Returns
        -------
        str
            Theseus output

        """
        output = subprocess.check_output(
            ["theseus", *pdbs_filename], stderr=subprocess.PIPE, universal_newlines=True
        )
        return output

    def _run_alignment(self):
        """
        Run MUSCLE

        Returns
        -------
        str
            Output of MUSCLE
        """
        output = subprocess.check_output(
            [
                self.alignment_app,
                "-maxiters",
                str(self.max_iterations),
                "-in",
                self.fastafile,
                "-out",
                self.alignment_file,
                "-clwstrict",
                "-verbose",
                "-log",
                "muscle.log",
            ],
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        _logger.info(output)
        return output

    def _run_theseus_different(self, pdbs_filename) -> str:
        """
        Superpose different sequences with Theseus based on the sequence alignment

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames

        Returns
        -------
        str
            Theseus output
        """
        self._get_fasta(pdbs_filename)
        self._concatenate_fasta(pdbs_filename)
        self._filemap(pdbs_filename)
        seq_alignment_output = self._run_alignment()

        self._parse_alignment(seq_alignment_output)
        _logger.info(seq_alignment_output)

        output = subprocess.check_output(
            ["theseus", "-f", "-M", self.filemap_file, "-A", self.alignment_file, *pdbs_filename],
            universal_newlines=True,
        )
        return output

    def _get_superposed_models(self, pdbs_filename) -> list:
        """
        Get the superposed models

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames

        Returns
        -------
        list
            list of superposed pdb models
        """
        superposed_pdb_filename = []
        for pdb in pdbs_filename:
            superposed_pdb_filename.append(f"theseus_{pdb}")

        self._strip_remark_lines(superposed_pdb_filename)

        superposed_pdb_models = [atomium.open(pdb_id).model for pdb_id in superposed_pdb_filename]
        return superposed_pdb_models

    def _strip_remark_lines(self, pdb_filenames) -> None:
        """
        Remove the "REMARK" lines in the pdb files

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames
        """
        for pdb in pdb_filenames:
            with open(pdb, "r") as infile:
                lines = infile.readlines()
            with open(pdb, "w") as outfile:
                for line in lines:
                    if not line.startswith("REMARK") or line.startswith("NUMMDL"):
                        outfile.write(line)

    def _get_transformation_matrix(self) -> dict:
        """
        Get the rotation matrix and translation vector

        Returns
        -------
        dict
            Rotation matrix and translation vector
        """
        translations, rotations = {}, {}
        with open("theseus_transf.txt", "r") as infile:
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
        return matrices

    def _parse_alignment(self, output) -> str:
        """
        Parse the output from the MSA program (muscle, by default)

        Parameter
        ---------
        output : bytes

        Returns
        -------
        str
            Output from MUSCLE
        -------
        """
        return output

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

    def run_theseus_different_no_superposition(self, structures) -> dict:
        """
        Calculate statistics (don't superposition)

        Parameter
        ---------
        pdbs_filename : list
            list of pdb filenames

        Returns
        -------
        dict
            As returned by ``._parse_superposition(output)``.
        """
        with enter_temp_directory(remove=False) as (cwd, tmpdir):
            _logger.info("All files are located in: %s", tmpdir)

            pdbs_filename = self._create_pdb(structures)
            self._get_fasta(pdbs_filename)
            self._concatenate_fasta(pdbs_filename)
            self._filemap(pdbs_filename)
            seq_alignment_output = self._run_alignment()

            self._parse_alignment(seq_alignment_output)
            _logger.info(seq_alignment_output)

            output = subprocess.check_output(
                [
                    "theseus",
                    "-I",
                    "-f",
                    "-M",
                    self.filemap_file,
                    "-A",
                    self.alignment_file,
                    *pdbs_filename,
                ],
                universal_newlines=True,
            )
        results = self._parse_superposition(output)
        return results
