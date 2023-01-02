from opencadd.chemical.composition import ProteinComposition
from opencadd.spacetime.volume import ToxelVolume

class BindingSite:
    def __init__(
            self,
            composition: ProteinComposition,
            volume: ToxelVolume,

    ):