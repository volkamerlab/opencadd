"""
Functions and routines for easier communication with PyMOL.
"""


# Standard library
from typing import Sequence, Tuple
from pathlib import Path
# 3rd party
from pymol import cmd
import numpy as np


def get_target_atoms_within_radius(
        filepath_target: Path,
        coords_center: Sequence[float, float, float],
        radius: float,
        properties: Sequence[str] = ("name", "resn"),
        reinitialize: bool = False,
):
    """
    Get the properties for all atoms in a target structure that are within a certain radius of a
    given point.

    Parameters
    ----------
    filepath_target : pathlib.Path
        Path to the target structure file.
    coords_center : Sequence[float, float, float]
        Coordinates (x, y, z) of the center of a sphere, to select all target atoms within.
    radius : float
        Radius of the sphere.
    properties : Sequence[str]
        A sequence of strings corresponding to the atom-property names defined in PyMol.
        For example, 'name' is the name of the atom, in PyMol notation, and `resn` is the name
        of the amino acid residue, to which that atom belongs.
    reinitialize : bool
        If True, PyMol will be reinitialized once before, and once after running this function,
        otherwise, no reinitialization will take place. Without reinitialization, everything
        done in PyMol during this process will remain as the active session, meaning future
        function calls to this or other PyMol functions will modify the same PyMol session,
        if those future calls also have this parameter set to False. Depending on your routine,
        both options can come in handy.

    Returns
    -------
    numpy.ndarray
        A 2-dimensional array of shape (n, k), where 'n' is equal to the number of found atoms,
        and 'k' is equal to the number of queried properties.
        Thus, the element at index (i, j) corresponds to the property at index 'j' in
        `properties`, for the 'i'th atom found by PyMol.
    """

    # TODO: This can most likely be done much faster if calculated internally using distance
    #  matrices between selected points and protein atoms.
    if reinitialize:
        cmd.reinitialize()
    # Load target structure from file.
    cmd.load(object="target_structure", filename=filepath_target)
    # Create pseudo-atom at the given center position.
    cmd.pseudoatom(object="center", pos=coords_center)
    # Select all atoms in target that are within the given radius to the center position.
    cmd.select(
        name="target_atoms_within_radius",
        selection=f"(center around {radius}) and target_structure",
    )
    # Iterate over selected atoms in PyMol, and extract values of given properties into a list.
    # TODO: define an internal convention to encode PyMol output strings into ints that can be
    #  optimally stored and accessed by numpy.
    selected_atoms_data = []
    cmd.iterate_state(
        state=1, # TODO: why 1? is it always 1?
        selection="target_atoms_within_radius",
        expression=f"selected_atoms_data.append([{', '.join(properties)}])",
        space={"selected_atoms_data": selected_atoms_data}
    )
    return np.array(selected_atoms_data)
