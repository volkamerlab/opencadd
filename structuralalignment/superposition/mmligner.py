class MMLignerWrapper(CommandLineWrapper):

    def __init__(self, option_a, option_b=None, ...):

        self.option_a = option_a
        self.option_b = option_b

    def _calculate(self, structures, **kwargs):
        # TODO: conda recipe for mmligner so executable is always there
        output = subprocess.check_output([
            self.executable,
            '--structure_a', structures[0].to_pdb(),
            '--structure_b', structures[1].to_pdb(),
            '--threshold', self.threshold,
            ...
        ])

        return self._parse(output)

    def _parse(self, output):
        return {
            'rmsd': ...,
            'score': ...,
            'metadata': {}
        }