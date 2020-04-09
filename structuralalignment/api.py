from .superposition.theseus import TheseusAligner
from .superposition.mmligner import MMLignerAligner
from .superposition.matchmaker import MatchMakerAligner

METHODS = {"theseus": TheseusAligner, "mmligner": MMLignerAligner, "matchmaker": MatchMakerAligner}


def align(*structures, method=TheseusAligner, **kwargs):
    """
    Main entry point for our project

    Parameters
    ----------
    structures : list of atomium.Model objects
        First one will be the targer to which the rest are aligned
    method : callable
        Usually a subclass of BaseAligner. This will be passed
        **kwargs

    Returns
    -------
    dict
        superposed models
        rmsd
        metadata
    """
    aligner = method(**kwargs)
    result = aligner.align(*structures)
    assert (
        len(result) == len(structures) - 1
    ), "Number of superpositions does not match number of structures"
    return result
