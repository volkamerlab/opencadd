import subprocess  # probably not supposed to be imported here?
from base import CommandLineWrapper


class MMLignerWrapper(CommandLineWrapper):
    """
    childclass of "CommandLineWrapper" that wrapps the mmligner to supperpose two protein structures"
    """

    def __init__(self):
        self.executable = "PATH"  # where is the executable located?

    def _calculate(self, structures, **kwargs):
        """
        calculates the superposition of two protein structures

        is called by CommandLineWrapper.calculate()

        Parameters
        ----------
        structures: [array like, array like]
            sequences of two protein structures of same length

        Returns
        -------
        self._parse(output)
            calls self._parse() to retrieve information about rmsd, score and metadata from output

        """
        # TODO: conda recipe for mmligner so executable is always there
        output = subprocess.check_output([
            self.executable,
            structures[0].to_pdb(),  # to_pdb() is provided by a different class?
            structures[1].to_pdb(),
            "--superpose"
        ])

        return self._parse(output)

    def _parse(self, output):
        """
        retrieves rmsd, score and metadata from the output of the mmligner subprocess

        Parameters
        ----------
        output: bytes
            string of bytes containing the stdout of the mmligener call

        Returns
        -------
        dict
            {
                "rmsd": float
                "score": float
                "metadata": ?
            }
            dictionary containing the rmsd value and the score of the superpoesed structures as well as the metadata
        """

        for line in output.splitlines():
            if line.startswith(b"RMSD"):  # maybe use "startswith" for speedup
                rmsd = float(line.split()[2])
            elif line.startswith(b"Coverage"):
                coverage = float(line.split()[2])

        return {
            'rmsd': rmsd,
            'score': coverage,  # score == coverage?
            'metadata': {}  # what's supposed to be returned as "metadata"?
        }
