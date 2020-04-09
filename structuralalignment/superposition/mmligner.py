import subprocess
import atomium
import biotite
import biotite.sequence.io.fasta as fasta

from base import CommandLineWrapper
from utils import enter_temp_directory


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

    Collier, J.H., Allison, L., Lesk A.M., Stuckey, P.J., Garcia de la Banda , M., Konagurthu, A.S. (2017)
    Statistical inference of protein structural alignments using information and compression.
    Bioinformatics, 33(7), 1005-1013

    Wallace,C.S. and Boulton,D.M. (1968) An information measure for classification.
    Comput. J., 11, 185–194.

    Wallace,C.S. (2005) Statistical and Inductive Inference Using MinimumMessage Length.
    Information Science and Statistics. SpringerVerlag, New York, NY.

    Shannon,C.E. (1948) A mathematical theory of communication.
    Bell Syst.Tech. J., 27, 379–423.
    """

    def __init__(self):
        # self.executable = "../../../../mmligner_1.0.2/bin/mmligner64.exe"
        self.executable = "C:/Universität/Softwarepraktikum/mmligner_1.0.2/bin/mmligner64.exe"

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

        with enter_temp_directory() as (cwd, tmpdir):

            path1, path2 = self._edit_pdb(structures)

            output = subprocess.check_output(
                [self.executable, path1, path2, "-o", "temp", "--superpose"]
            )

            result = self._parse(output)

        return result

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
                ivalue = float(line.split()[4])

        alignment = fasta.FastaFile()

        return {
            "rmsd": rmsd,
            "score": ivalue,
            "metadata": {"alignment": alignment.read("./temp__1.afasta")},
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
        dictionary containing the rmsd value and the score of the alignment for the given two structures as well as the
        metadata by calling self._parse(output)
        """
        assert len(structures) == 2

        output = subprocess.check_output(
            [
                self.executable,
                structures[0].to_pdb(),
                structures[1].to_pdb(),
                "--ivalue",
                alignment.to_fasta(),
            ]
        )

        return self._parse(output)

    def _edit_pdb(self, structures, path=["./structure1.pdb", "./structure2.pdb"]):
        """
        function to write atomium protein models to pdb readable by mmligner

        Parameters
        ----------
        structures: [array like, array like]
            two protein structures

        path: [str, str], Optional=["structure1.pdb, "structure2.pdb"]
            Path where the pdbs should be written

        Returns
        -------
        str, str
            Paths of both structures
        """
        assert len(path) == 2

        structures[0].save(path[0])
        structures[1].save(path[1])

        for i in range(len(path)):
            pdb = []
            with open(path[i], "r") as infile:
                assert infile.mode == "r"

                pdb = infile.readlines()
                for j in range(1, len(pdb)):
                    if pdb[j].startswith("TER"):
                        pdb[j] = pdb[j].split()[0] + "    " + pdb[j - 1].split()[1] + "\n"

            self._write_pdb(path[i], pdb)

        return path[0], path[1]

    def _write_pdb(self, path, pdb):
        """
        function to write atomium protein models to pdb readable by mmligner

        Parameters
        ----------
        path: str
            Path where the pdb should be written

        pdb: array like
            edited pdb file

        Returns
        -------
        none
        """
        with open(path, "w") as outfile:
            assert outfile.mode == "w"

            for line in pdb:
                outfile.write(line)

        return


if __name__ == "__main__":
    """
    Script for debugging.
    """
    structures = []
    structures.append(atomium.fetch("4u3y").model)
    structures.append(atomium.fetch("4u40").model)

    Aligner = MMLignerWrapper()

    # Aligner._edit_pdb(structures)

    result = Aligner.calculate(structures)
    print(f"Alignment has a RMSD of: {result['rmsd']} and a Score (IValue) of: {result['score']} ")
    # print("Alignment data:")
    # print(result["metadata"]["alignment"])
    # for header, string in result["metadata"]["alignment"].items():
    #    print("Header:", header)
    #    print(len(string))
    #    print("Sequence:", string[:50], "...")
    #    print("Sequence length:", len(string))
