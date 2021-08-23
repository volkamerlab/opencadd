"""
Defines easy programmatic access for any entry point
"""

from .engines.theseus import TheseusAligner
from .engines.mmligner import MMLignerAligner
from .engines.mda import MDAnalysisAligner
from ..core import Structure


METHODS = {
    "theseus": TheseusAligner,
    "mmligner": MMLignerAligner,
    "mda": MDAnalysisAligner,
}


def align(structures, method=TheseusAligner, **kwargs):
    """
    Main entry point for our project

    Parameters
    ----------
    structures : list of opencadd.core.Structure objects
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
    assert all(isinstance(s, Structure) for s in structures)
    reference, *mobiles = structures
    results = []
    for mobile in mobiles:
        # only take the first models of the pdb files, this ensures that mda and theseus are working consitently
        # comparing all models could provide better results, but would be very inefficient
        # (e.g. 25 models would mean 25 times the computing time)
        result = aligner.calculate([reference.models[0], mobile.models[0]])
        results.append(result)

    return results
