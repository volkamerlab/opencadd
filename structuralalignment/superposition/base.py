class CommandLineWrapper:
    def calculate(self, structures):
        assert len(structures) == 2
        return self._calculate(structures)

    def _calculate(self, structures, **kwargs):
        raise NotImplementedError("Reimplement in your subclass")


class BaseAligner:
    def align(self, *args, **kwargs):
        print("Hey, it's working. I am supposed to do something useful via a subclass.")

    def compute_sequence_alignment(self, algorithm, **kwargs):
        """
        Example implementation of first step
        """
        sequence_1 = self.structures[0].sequence()
        sequence_2 = self.structures[1].sequence()
        sequence_aligner = algorithm(**kwargs)
        alignment = sequence_aligner(sequence_1, sequence_2, **kwargs)
        return alignment

    def compute_structural_overlap(self, alignment, algorithm, **kwargs):
        structural_aligner = algorithm(**kwargs)
        rmsd_before = self.compute_rmsd(self.structures[0], self.structures[1])
        overlapped_1, overlapped_2 = structural_aligner(
            self.structures[0], self.structures[1], alignment
        )
        rmsd_after = self.compute_rmsd(overlapped_1, overlapped_2)
        logger.log("RMSD before and after: %d, %d", rmsd_before, rmsd_after)
        return overlapped_1, overlapped_2, rmsd_after
