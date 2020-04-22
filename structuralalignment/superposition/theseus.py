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

.. todo::

    - Add literature references to this docstring.
    - Maybe parse more information from output files
    - Try different alignment tool (biotite)
    - Improve documentation
    - Delete debug code
"""

import subprocess
from pathlib import Path
import logging
import atomium

import numpy as np

from structuralalignment.superposition.base import BaseAligner

# from ..sequences import get_alignment_fasta
import biotite.sequence.io.fasta as fasta
from structuralalignment.utils import enter_temp_directory


_logger = logging.getLogger(__name__)


class TheseusAligner(BaseAligner):

    """
    Superpose structures with different sequences but similar structures

    Parameters
    ----------
    max_iterations : int
        number of iterations for muscle

    Examples
    --------

    .. todo::

        Move these two examples to ``docs/examples`` or ``docs/tutorials``

    * Multiple structures of the cytochrome c protein from different species
    * Multiple mutated structures of hen egg white lysozyme.
    """

    def __init__(self, max_iterations: int = 32):
        self.max_iterations = max_iterations

        self.fastafile = "theseus.fasta"
        self.filemap_file = "theseus.filemap"
        self.alignment_file = "theseus.aln"
        self.alignment_file_biotite = "theseus_biotite.aln"
        self.alignment_app = "muscle"

    def _create_pdb(self, structures):
        """
        Create files with atomium models

        Parameter
        ---------
        structures : list
            list of atomium models
        """
        fetched_pdbs_filename = []
        for index, pdbs in enumerate(structures):
            pdb_filename = f"{index}.pdb"
            pdbs.save(pdb_filename)
            fetched_pdbs_filename.append(pdb_filename)
        return fetched_pdbs_filename

    def _get_fasta(self, pdbs_filename):
        """
        Use Theseus to create fasta files

        .. todo::
            We can probably generate this from atomium directly
        """
        return subprocess.check_output(["theseus", "-f", "-F", *pdbs_filename])

    def _concatenate_fasta(self, pdbs_filename):
        """
        Concatenate the fasta files created with Theseus into one multiple sequence fasta file
        """
        with open(self.fastafile, "w") as outfile:
            for fname in pdbs_filename:
                with open(f"{fname}.fst") as infile:
                    for line in infile:
                        outfile.write(line)

    def _filemap(self, pdbs_filename):
        """
        Create filemap for Theseus
        """
        with open(self.filemap_file, "w") as outfile:
            for fname in pdbs_filename:
                outfile.write(f"{fname} {fname}\n")

    # pylint: disable=arguments-differ
    def _calculate(self, structures, identical: bool = False, **kwargs):
        """
        Align the sequences with an alignment tool (``muscle``)

        .. todo::

            Can we use a different alignment tool, preferrably in Python? E.g. biotite?


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

            - ``rmsd`` (float): RMSD Value of the alignment
            - ``score`` (float)
                    ivalue of the alignment. The smaller the better
            - ``metadata`` (?)
        """

        with enter_temp_directory(remove=False) as (cwd, tmpdir):
            _logger.info("DEBUG: Running in %s -- TODO: DELETE BEFORE MERGE --", tmpdir)

            pdbs_filename = self._create_pdb(structures)

            if identical:
                superposition_output = self._run_theseus_identical(pdbs_filename)
            else:
                superposition_output = self._run_theseus_different(pdbs_filename)
                _logger.info(superposition_output)
                _logger.info("-- OUTPUT WE NEED TO PARSE BEFORE DELETING TEMPFILES --")
                _logger.info("\n".join([str(item) for item in Path(tmpdir).glob("theseus_*")]))
            superposed_pdb_models = self._get_superposed_models(pdbs_filename)
            transformation_matrix = self._get_transformation_matrix()
        results = self._parse_superposition(superposition_output)
        results["superposed"] = superposed_pdb_models
        results["metadata"]["transformation"] = transformation_matrix
        return results

    def _run_theseus_identical(self, pdbs_filename):
        """
        Superpose identical sequences with Theseus
        """
        output = subprocess.check_output(
            ["theseus", *pdbs_filename], stderr=subprocess.PIPE, universal_newlines=True
        )
        return output

    def _run_alignment(self):
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

    def _run_alignment_biotite(self):
        # TODO: seq_strings = list(fasta_file.values()) // AttributeError: 'str' object has no attribute 'values'
        fasta_file = fasta.FastaFile()
        fasta_file.read(self.fastafile)
        alignment = fasta.get_alignment(fasta_file, additional_gap_chars=("-",))
        with open(self.alignment_file_biotite, "w") as outfile:
            outfile.write(alignment)
        return alignment

    def _run_theseus_different(self, pdbs_filename):
        """
        Superpose different sequences with Theseus based on the sequence alignment
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

    def _get_superposed_models(self, pdbs_filename):
        """
        .. todo:

            Get the coordinates of superposed_pdb_models and add them to the original
            Atomium Models or maybe a copy of those.

        """
        superposed_pdb_filename = []
        for pdb in pdbs_filename:
            superposed_pdb_filename.append(f"theseus_{pdb}")

        self._strip_remark_lines(superposed_pdb_filename)

        superposed_pdb_models = [atomium.open(pdb_id).model for pdb_id in superposed_pdb_filename]
        return superposed_pdb_models

    def _strip_remark_lines(self, pdb_filenames):
        for pdb in pdb_filenames:
            with open(pdb, "r") as infile:
                lines = infile.readlines()
            with open(pdb, "w") as outfile:
                for line in lines:
                    if not line.startswith("REMARK") or line.startswith("NUMMDL"):
                        outfile.write(line)

    def _get_transformation_matrix(self):
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

    def _parse_alignment(self, output):
        """
        Parse the output from the MSA program (muscle, by default)

        Parameters
        ----------
        output : bytes
        """
        return output

    def _parse_superposition(self, output):
        """
        Parse the output from theseus itself

        Parameters
        ----------
        output : bytes
        """
        print(output)
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
                # TODO: check for 1.13e-02 value
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
            },  # TODO: add info from top
            # TODO:add residues
        }


# TODO: REMOVE
if __name__ == "__main__":
    pdb_ids = ["6HG4", "6HG9"]

    models = [atomium.fetch(pdb_id).model for pdb_id in pdb_ids]

    theseus = TheseusAligner()

    theseus._calculate(models, False)
