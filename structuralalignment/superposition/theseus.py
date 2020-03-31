import subprocess
from pathlib import Path

import atomium

from structuralalignment.superposition.base import CommandLineWrapper
from structuralalignment.utils import enter_temp_directory


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
    - Can we unify both classes in a single one? Mode of operation
      would be configured with a keyword (e.g. `identical=True`).
"""


class Theseus(CommandLineWrapper):
    """
    A mode for superpositioning macromolecules with
    identical sequences and numbers of residues,
    for instance, multiple models in an NMR family
    or multiple structures from different crystal
    forms of the same protein.
    """

    def _calculate(self, structures, **kwargs):
        """
        Superpose structures with identical sequences with the tool Theseus

        Parameter
        ---------
        structures : list
            list of atomium objects

        **kwargs : dictionary
            optional parameter
        """
        fetched_pdbs = []
        fetched_pdbs_filename = []

        for pdb in structures:
            fetched_pdbs.append(pdb)
        with enter_temp_directory():
            for pdb in fetched_pdbs:
                pdb_filename = pdb.code + "." + pdb.filetype
                pdb.model.save(pdb_filename)
                fetched_pdbs_filename.append(pdb_filename)

            output = subprocess.check_output(["theseus"] + fetched_pdbs_filename)

        return self._parse(output)

    def _parse(self, output):
        """
        Parse the output
        """
        return output


class TheseusAlign(CommandLineWrapper):

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
        Create files with atomium objects

        Parameter
        ---------
        structures : list
            list of atomium objects

        .. todo::
            This method should be private and not set state,
            just return temporary filenames
        """
        self.fetched_pdbs = structures
        self.fetched_pdbs_filename = []

        for pdbs in self.fetched_pdbs:
            pdb_filename = f"{pdbs.code}.pdb"
            pdbs.model.save(pdb_filename)
            self.fetched_pdbs_filename.append(pdb_filename)

    def _get_fasta(self):
        """
        Use Theseus to create fasta files

        .. todo::
            We can probably generate this from atomium directly
        """
        return subprocess.check_output(["theseus", "-f", "-F", *self.fetched_pdbs_filename])

    def _concatenate_fasta(self):
        """
        Concatenate the fasta files created with Theseus into one multiple sequence fasta file
        """
        filenames = self.fetched_pdbs_filename
        with open(self.fastafile, "w") as outfile:
            for fname in filenames:
                with open(f"{fname}.fst") as infile:
                    for line in infile:
                        outfile.write(line)

    def _filemap(self):
        """
        Create filemap for Theseus
        """
        filenames = self.fetched_pdbs_filename
        with open(self.filemap_file, "w") as outfile:
            for fname in filenames:
                outfile.write(f"{fname} {fname}\n")

    def _calculate(self, structures, **kwargs):
        """
        Align the sequences with an alignment tool (``muscle``)

        .. todo::

            Can we use a different alignment tool, preferrably in Python? E.g. biotite?

        Parameters
        ----------
        structures : list
            list of atomium objects
        """
        with enter_temp_directory(remove=False) as (cwd, tmpdir):
            print("DEBUG: Running in", tmpdir, " -- TODO: DELETE BEFORE MERGE --")
            self._create_pdb(structures)
            self._get_fasta()
            self._concatenate_fasta()
            self._filemap()
            seq_alignment_output = subprocess.check_output(
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
                universal_newlines=True,
            )
            self._parse_alignment(seq_alignment_output)
            print(seq_alignment_output)
            superposition_output = self._run_theseus()
            print(superposition_output)
            print("-- OUTPUT WE NEED TO PARSE BEFORE DELETING TEMPFILES --")
            print(*list(Path(tmpdir).glob("theseus_*")), sep="\n")
        return self._parse_superposition(superposition_output)

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
        return output

    def _run_theseus(self):
        """
        Superpose with Theseus based on the sequence alignment
        """
        output = subprocess.check_output(
            [
                "theseus",
                "-f",
                "-M",
                self.filemap_file,
                "-A",
                self.alignment_file,
                *self.fetched_pdbs_filename,
            ],
            universal_newlines=True,
        )
        return output


if __name__ == "__main__":
    import sys
    import atomium

    pdb_ids = sys.argv[1:]
    models = [atomium.fetch(pdb_id) for pdb_id in pdb_ids]
    theseus = TheseusAlign()
    theseus.calculate(models)
