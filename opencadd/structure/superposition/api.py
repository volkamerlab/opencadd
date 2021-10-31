# was zum teufel habe ich da gemacht

"""
Defines easy programmatic access for any entry point
"""

from MDAnalysis import _SELECTION_WRITERS
from .engines.theseus import TheseusAligner
from .engines.mmligner import MMLignerAligner
from .engines.mda import MDAnalysisAligner
from ..core import Structure


METHODS = {
    "theseus": TheseusAligner,
    "mmligner": MMLignerAligner,
    "mda": MDAnalysisAligner,
}

def align(structures, user_select, method=TheseusAligner, only_backbone=False,**kwargs):
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


    if user_select:
        # only take the first models of the pdb files, this ensures that mda and theseus are working consitently
        # comparing all models could provide better results, but would be very inefficient
        # (e.g. 25 models would mean 25 times the computing time)
            # selection always returns an AtomGroup. Alternative, which is not recommended, 
            # because the aligners can handle AtomGroups aswell: 
            # Structure.from_atomgroup(reference.models[0].select_atoms(f"{user_select}").residues.atoms)
        reference = reference.models[0]
        reference_selection = reference.select_atoms(f"{user_select[0]}")
        for mobile in mobiles:
            mobile = mobile.models[0]
            mobile_selection = mobile.select_atoms(f"{user_select[1]}")
            result = aligner.calculate(structures = [reference, mobile], selections = [reference_selection, mobile_selection])
            # size of selections
            reference_size = len(reference_selection.residues)
            mobile_size = len(mobile_selection.residues)
            result["metadata"]["reference_size"] = reference_size
            result["metadata"]["mobile_size"] = mobile_size
            results.append(result)

    else:
        reference = reference.models[0]
        for mobile in mobiles:
            mobile = mobile.models[0]
            result = aligner.calculate(structures = [reference, mobile], selections = [])
            reference_size = len(reference.residues)
            mobile_size = len(mobile.residues)
            result["metadata"]["reference_size"] = reference_size
            result["metadata"]["mobile_size"] = mobile_size
            results.append(result)

    return results
