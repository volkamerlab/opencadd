import atomium
import subprocess
import logging
logger = logging.getLogger(__name__)

# usage: Theseus(["3dk3"])
class Theseus():
    """
    Superpose structures with identical sequences with the tool Theseus
    """

    def align(*structures):
        """
        Fetch structures with atomium
        Save the fetched PDB files and run the theseus command with subprocess

        Parameter
        ---------
        structures : list
            list of PDB file names
        """
        fetched_pdbs = []
        fetched_pdbs_filename = []

        for pdb in structures:
            fetched_pdbs.append(atomium.fetch(pdb))

        for pdbs in fetched_pdbs:
            pdb_filename = pdbs.code + "." + pdbs.filetype
            pdbs.model.save(
                pdb_filename
            )
            fetched_pdbs_filename.append(
                pdb_filename
            )

        output = subprocess.check_output(["theseus"] + fetched_pdbs_filename)
        return _log(output)

    def _log(self, output):
        """
        Log the output
        """
        logger.info(output)

"""
usage: call class first
a=TheseusAlign()
a.align("3DK3")
"""
class TheseusAlign:
    """
    Superpose structures with different sequences but similar structures
    """

    def __init__(self):
        self.fastafile = "theseus.fasta"
        self.filemap_file = "theseus.filemap"
        self.alignment_file = "theseus.aln"
        self.alignment_app = "muscle"

    def fetch_pdb(self, structures):
        """
        Fetch structures with atomium

        Parameter
        ---------
        structures : list
            list of PDB file names
        """
        self.fetched_pdbs = []
        self.fetched_pdbs_filename = []

        for pdb in structures:
            self.fetched_pdbs.append(atomium.fetch(pdb))

        for pdbs in self.fetched_pdbs:
            pdb_filename = pdbs.code + "." + pdbs.filetype
            pdbs.model.save(
                pdb_filename
            )
            self.fetched_pdbs_filename.append(
                pdb_filename
            )

    def get_fasta(self):
        """
        Use Theseus to make fasta files
        """

        output = subprocess.check_output(
            ["theseus", "-f", "-F"] + self.fetched_pdbs_filename
        )
        self._log(output)

    def concatenate_fasta(self):
        """
        Concatenate the fasta files made with Theseus into one multiple sequence fasta file
        """
        filenames = self.fetched_pdbs_filename
        with open(self.fastafile, "w") as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def filemap(self):
        """
        Create filemap for Theseus
        """
        filenames = self.fetched_pdbs_filename
        with open(self.filemap_file, "w") as outfile:
            for fname in filenames:
                outfile.write(fname + " " + fname)

    def align(self, *structures):
        """
        Align the sequences with an alignment tool
        """
        self.fetch_pdb(structures)
        self.get_fasta()
        self.concatenate_fasta()
        self.filemap()
        output = subprocess.check_output(
            [
                self.alignment_app,
                "-maxiters",
                "32",
                "-in",
                self.fastafile,
                "-out",
                self.alignment_file,
                "-clwstrict",
            ]
        )
        self.execute()
        return self._log(output)

    def _log(self, output):
        """
        Parse the output
        """
        logger.info(output)

    def execute(self):
        """
        Superpose with Theseus based on the sequence alignment
        """
        output = subprocess.check_output(
            [
                "-f",
                "-M",
                self.filemap_file,
                "-A",
                self.alignment_file,
                self.fetched_pdbs_filename,
            ]
        )
        self._log(output)

