"""
Defines easy programmatic access for any entry point
"""

from MDAnalysis.exceptions import NoDataError


from .engines.theseus import TheseusAligner
from .engines.mmligner import MMLignerAligner
from .engines.mda import MDAnalysisAligner
from ..core import Structure


METHODS = {
    "theseus": TheseusAligner,
    "mmligner": MMLignerAligner,
    "mda": MDAnalysisAligner,
}


def align(structures, user_select=False, method=TheseusAligner, **kwargs):
    """
    Main entry point for our project

    Parameters
    ----------
    structures : list of opencadd.core.Structure objects
        First one will be the target to which the rest of the structures are aligned.
    user_select: list of MDAnalysis selection strings
        Provided by user in the CLI (default: False).
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

    def _universe_to_atomgroup(universe):
        """
        Universes can come with or without models. For downstream processing, cast all universes 
        to AtomGroups.

        Parameters
        ----------
        universe : MDAnalysis.core.universe.Universe
            Input structure with or without models.
        
        Returns
        -------
        MDAnalysis.core.groups.AtomGroup
            Structure as AtomGroup (keep only first model if multiple models available).
        """
        try:
            atom_group = universe.models[0]
        except NoDataError:
            atom_group = universe.atoms
        return atom_group

    # only take the first models of the pdb files, this ensures that mda and theseus are working consistently
    # comparing all models could provide better results, but would be very inefficient
    # (e.g. 25 models would mean 25 times the computing time)
    if user_select:
        # selection always returns an AtomGroup. Alternative, which is not recommended,
        # because the aligners can handle AtomGroups aswell:
        # Structure.from_atomgroup(reference.models[0].select_atoms(f"{user_select}").residues.atoms)
        reference = _universe_to_atomgroup(reference)
        reference_selection = reference.select_atoms(f"{user_select[0]}")
        for mobile in mobiles:
            mobile = _universe_to_atomgroup(mobile)
            mobile_selection = mobile.select_atoms(f"{user_select[1]}")
            result = aligner.calculate(
                structures=[reference, mobile], selections=[reference_selection, mobile_selection]
            )
            # size of selections
            reference_size = len(reference_selection.residues)
            mobile_size = len(mobile_selection.residues)
            result["metadata"]["reference_size"] = reference_size
            result["metadata"]["mobile_size"] = mobile_size
            results.append(result)

    else:
        reference = _universe_to_atomgroup(reference)
        for mobile in mobiles:
            mobile = _universe_to_atomgroup(mobile)
            result = aligner.calculate(structures=[reference, mobile], selections=[])
            # size of whole structures
            reference_size = len(reference.residues)
            mobile_size = len(mobile.residues)
            result["metadata"]["reference_size"] = reference_size
            result["metadata"]["mobile_size"] = mobile_size
            results.append(result)

    return results
