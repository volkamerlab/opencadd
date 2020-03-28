import subprocess
from .base import CommandLineWrapper


class MMLignerWrapper(CommandLineWrapper):
    """
    subclass of "CommandLineWrapper" that wrapps the mmligner to superpose two protein structures

    The MMLigner (Collier et al., 2017) works by minimizing the ivalue of the alignment. The ivalue is based on
    the Minimum Message Length framework (Wallace and Boulton, 1968; Wallace, 2005), a Bayesian framework for
    statistical inductive inference. The ivalue represents the hypothetical minimum message length needed to transmit
    the computed alignment losslessly (Shannon, 1948).
    Using the ivalue measure, the algorithm creates crude-but-effective strucural alignments rapidly to act as seeds.
    These seeds are iteratively refined over an Expectation-Maximization approach using the I-value criterion.
    By substracting the ivalue from the null model, the statistical significance of the alignment can be computed. If the
    difference is positive, the alignment is significant.
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
                    RMSD Value of the alignment
                "score": float
                    ivalue of the alignment. The smaller the better
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
                    RMSD Value of the alignment
                "score": float
                    ivalue of the alignment. The smaller the better
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
                    RMSD Value of the alignment
                "score": float
                    ivalue of the alignment. The smaller the better
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
