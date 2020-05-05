"""
Defines easy programmatic access for any entry point
"""

from .superposition.theseus import TheseusAligner
from .superposition.mmligner import MMLignerAligner
from .superposition.matchmaker import MatchMakerAligner
from .core import Structure


METHODS = {
    "theseus": TheseusAligner,
    "mmligner": MMLignerAligner,
    "matchmaker": MatchMakerAligner,
}


def align(structures, method=TheseusAligner, **kwargs):
    """
    Main entry point for our project

    Parameters
    ----------
    structures : list of superposer.core.Structure objects
        First one will be the targer to which the rest are aligned
    method : BaseAligner-like
        Usually a subclass of BaseAligner. This will be passed ``**kwargs``. This class
        MUST define `.calculate()`.

    Returns
    -------
    dict
        superposed models
        rmsd
        metadata
    """
    aligner = method(**kwargs)
    reference, *mobiles = structures
    results = []
    for mobile in mobiles:
        result = aligner.calculate([reference, mobile])
        results.append(result)

    return results
