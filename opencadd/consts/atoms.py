from typing import Literal

import numpy as np

from opencadd._typing import ArrayLike


def van_der_waals_radius(
        atom_ids: ArrayLike,
        id_type: Literal["internal", "element_symbol"] = "internal"
) -> np.ndarray:
    """
    Get the Van der Waals radii of a set of given atoms.

    Parameters
    ----------
    atom_ids : opencadd.typing.ArrayLike
        1D array of atom identifiers.
    id_type : Literal["internal", "element_symbol"], optional, default: "internal"
        Type of atom identifiers.

    Returns
    -------
    numpy.ndarray
        1D array of van der Waals radii corresponding to the input array.

    References
    ----------
    https://github.com/rdkit/rdkit/issues/1831

    """
    raise NotImplementedError





