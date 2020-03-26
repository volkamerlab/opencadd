import subprocess
from .base import CommandLineWrapper


class MMLignerWrapper(CommandLineWrapper):
    """
    subclass of "CommandLineWrapper" that wrapps the mmligner to superpose two protein structures
    """

    def __init__(self):
        self.executable = mmligner

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
        dict
            {
                "rmsd": float
                "score": float
                "metadata": ?
            }
            dictionary containing the rmsd value and the score of the superpoesed structures as well as the metadata. Parsed by self._parse(output).

        """
        # TODO: conda recipe for mmligner so executable is always there
        output = subprocess.check_output([
            self.executable,
            structures[0].to_pdb(),
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
            if line.startswith(b"RMSD"):
                rmsd = float(line.split()[2])

            # coverage may still be relevant?
            # elif line.startswith(b"Coverage"):
            #   coverage = float(line.split()[2])

            elif line.startswith(b"I(A & <S,T>)"):
                ivalue = float(line.split()[2])

        return {
            'rmsd': rmsd,
            'score': ivalue,
            'metadata': {}  # what's supposed to be returned as "metadata"?
        }

    def ivalue(self, structures, alignment):
        """
        computes the score and rmsd for a given alignment of two structures by calling mmligner as a subprocess

        Parameters
        ----------
        structures: [array like, array like]
            sequences of two protein structures of same length

        alignment: array like
            alignment of the given two sequences

        Returns
        -------
        dict
            {
                "rmsd": float
                "score": float
                "metadata": ?
            }
        dictionary containing the rmsd value and the score of the alignment for the given two structures as well as the metadata by calling self._parse(output)
        """
        assert len(structures) == 2

        output = subprocess.check_output([
            self.executable,
            structures[0].to_pdb(),
            structures[1].to_pdb(),
            "--ivalue", alignment.to_fasta()
        ])

        return self._parse(output)
