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
        return self._parse_superposition(superposition_output)

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
            ],
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        return output

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
        for line in output.splitlines():
            if "Classical" in line:
                rmsd = float(line.split()[5])

        return {
            "superposed": None,  # TODO: Add the superposed models here!
            "scores": {"rmsd": rmsd},  # TODO: More scores from output?
            "metadata": {},  # TODO: See what interesting extra info we have in the output
        }
