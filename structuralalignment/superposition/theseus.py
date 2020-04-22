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

from structuralalignment.superposition.base import BaseAligner
#from ..sequences import get_alignment_fasta
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
        return self._parse_superposition(superposition_output, superposed_pdb_models, transformation_matrix)

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
        read_fasta = fasta_file.read(self.fastafile)
        alignment = fasta.get_alignment(read_fasta, additional_gap_chars = ('-',))
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
                    if not line.startswith('REMARK') or line.startswith('NUMMDL'):
                        outfile.write(line)

    def _get_transformation_matrix(self):
        matrix = {}
        with open("theseus_transf.txt", "r") as infile:
            lines = infile.readlines()
            for line in lines:
                if "R:" in line:
                    model = line.split(":")[0]
                    matrices = line.split(":")[1]
                    matrix[model] = matrices
        return matrix

    def _parse_alignment(self, output):
        """
        Parse the output from the MSA program (muscle, by default)

        Parameters
        ----------
        output : bytes
        """
        return output

    def _parse_superposition(self, output, superposed_pdb_models, transformation_matrix):
        """
        Parse the output from theseus itself

        Parameters
        ----------
        output : bytes
        """
        for line in output.splitlines():
            if "Classical" in line:
                rmsd = float(line.split()[5])
            if "Least-squares" in line:
                least_squares = float(line.split()[3])
            if "Maximum" in line:
                maximum_likelihood = float(line.split()[4])
            if "Marginal" in line:
                log_marginal_likelihood = float(line.split()[4])
            if "AIC" in line:
                aic = float(line.split()[2])
            if "BIC" in line:
                bic = float(line.split()[2])
            if "Omnibus" in line:
                omnibus_chi_square = float(line.split()[3])
            if "Hierarchical var" in line:
                hierarchical_var_chi_square = float(line.split()[5])
            if "Rotational" in line:
                rotational_translational_covar_chi_square = float(line.split()[5])
            if "Hierarchical minimum" in line:
                hierarchical_minimum_var_sigma = float(line.split()[5]) # TODO: check for 1.13e-02 value
            if "skewness" in line:
                skewness = float(line.split()[2])
            if "skewness Z-value" in line:
                skewness_z = float(line.split()[3])
            if "kurtosis" in line:
                kurtosis = float(line.split()[2])
            if "kurtosis Z-value" in line:
                kurtosis_z = float(line.split()[3])
            if "data pts" in line:
                data_pts = float(line.split()[4])
                free_params = float(line.split()[8])
                d_p = float(line.split()[11])
            if "Median" in line:
                median_structure = float(line.split()[4])
            if "N(total)" in line:
                n_total = float(line.split()[3])
                n_atoms = float(line.split()[6])
                n_structures = float(line.split()[9])
            if "Total rounds" in line:
                total_rounds = float(line.split()[3])
        return {
            "superposed": superposed_pdb_models, #TODO is this correct or should the models be in {}?
            "scores": {"rmsd": rmsd},
            "metadata": {"transformation": transformation_matrix},  # TODO: add info from top
            #TODO:add residues
        }


# TODO: REMOVE
if __name__ == "__main__":
    pdb_ids = ["6HG4", "6HG9"]

    models = [atomium.fetch(pdb_id).model for pdb_id in pdb_ids]

    theseus = TheseusAligner()

    theseus._calculate(models, False)
