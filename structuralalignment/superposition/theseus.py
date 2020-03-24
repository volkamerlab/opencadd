import atomium
import subprocess

# usage: Theseus(["3dk3"])
class Theseus:
    """
    Superpose structures with identical sequences with the tool Theseus
    """

    def __init__(self, pdbs):
        """
        Initialize variables and fetch PDB files with atomium

        Parameters
        ----------
        pdbs : array-like
            array of PDB file names
        """
        self.pdbs = pdbs
        self.fetched_pdbs = []
        for pdb in pdbs:
            self.fetched_pdbs.append(atomium.fetch(pdb))
        self.align()

    def align(self):
        """
        Save the fetched PDB files and run the theseus command with subprocess
        """

        self.fetched_pdbs_filename = []
        """
        for m in range(len(self.fetched_pdbs)):
            self.fetched_pdbs[m].model.save(self.pdbs[m])
        """
        for m in range(len(self.fetched_pdbs)):
            self.fetched_pdbs[m].model.save(
                self.pdbs[m] + "." + self.fetched_pdbs[m].filetype
            )
            self.fetched_pdbs_filename.append(
                self.pdbs[m] + "." + self.fetched_pdbs[m].filetype
            )

        output = subprocess.check_output(["theseus",] + self.fetched_pdbs_filename)
        return self._parse(output)

    def _parse(self, output):
        """
        Parse the output
        """
        print(output)


# usage: TheseusAlign(["3DK3"])
class TheseusAlign:
    """
    Superpose structures with different sequences but similar structures
    """

    def __init__(self, pdbs):
        self.pdbs = pdbs
        self.fetched_pdbs = []
        self.fetched_pdbs_filename = []

        self.fastafile = "theseus.fasta"
        self.filemap_file = "theseus.filemap"
        self.alignment_file = "theseus.aln"
        self.alignment_app = "/usr/local/bin/muscle"

        for pdb in pdbs:
            self.fetched_pdbs.append(atomium.fetch(pdb))

        self.get_fasta()
        self.concatenate_fasta()
        self.align()
        self.use_theseus()

    def get_fasta(self):
        """
        Use Theseus to make fasta files
        """
        for m in range(len(self.fetched_pdbs)):
            self.fetched_pdbs[m].model.save(
                self.pdbs[m] + "." + self.fetched_pdbs[m].filetype
            )
            self.fetched_pdbs_filename.append(
                self.pdbs[m] + "." + self.fetched_pdbs[m].filetype
            )

        output = subprocess.check_output(
            ["theseus", "-f", "-F"] + self.fetched_pdbs_filename
        )

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

    def align(self):
        """
        Align the sequences with an alignment tool
        """
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
        return self._parse(output)

    def _parse(self, output):
        """
        Parse the output
        """
        print(output)

    def use_theseus(self):
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
        self._parse(output)

