
# Standard library
from typing import Sequence, Union, Optional, Tuple, Literal
import subprocess
from pathlib import Path
# 3rd-party
import numpy as np
# Self
from opencadd._typing import PathLike, ArrayLike


def extract_field_values(
        map_filepath: PathLike,
        data_type: np.dtype = np.single,
) -> np.ndarray:
    """
    Extract the calculated grid point energies from a MAP output file.

    Parameters
    ----------
    map_filepath : pathlib.Path
        Filepath of a MAP file outputted by AutoGrid.
    data_type : numpy.dtype, Optional, default: numpy.single
        Numpy datatype of the output array. Default is 32-bit float (numpy.single).

    Returns
    -------
    grid : numpy.ndarray
        A 1-dimensional array containing the calculated energy values for each grid point.
        The grid points are ordered according to the nested loops z(y(x)), so the x-coordinate
        is changing fastest. The coordinate system is right-handed.

    Notes
    -----
    In a MAP file, the first 6 lines are headers. The energy values start at line 7,
    and are written one per line, until the end of file.

    Example of first few lines of a MAP file:
    ```
    GRID_PARAMETER_FILE vac1.nbc.gpf
    GRID_DATA_FILE 4phv.nbc_maps.fld
    MACROMOLECULE 4phv.new.pdbq
    SPACING 0.375
    NELEMENTS 50 50 80
    CENTER -0.026 4.353 -0.038
    125.095596
    123.634560
    116.724602
    108.233879
    ```

    The header `NELEMENTS` is the same as the input parameter `npts`, defined in function
    `create_gpf`.
    """
    with open(map_filepath, "r") as f:
        lines = f.readlines()
    return np.array(lines[6:]).astype(data_type)


def extract_grid_params(
        map_filepath: PathLike
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], float]:
    """
    Extract the AutoGrid parameters `gridcenter`, `npts` and `spacing` from a MAP file.

    Parameters
    ----------
    map_filepath : pathlib.Path
        Path to one of the calculated MAP files.

    Returns
    -------
    gridcenter, npts, spacing : tuple[tuple[float, float, float], tuple[float, float, float], float]

    See Also
    --------
    For a sample of MAP file contents, see function `create_grid_tensor_from_map_file`.
    """
    with open(map_filepath, "r") as f:
        lines = f.readlines()
    grid_spacing = float(lines[3].split()[1])
    npts = tuple(map(float, lines[4].split()[1:4]))
    grid_center = tuple(map(float, lines[5].split()[1:4]))
    return grid_center, npts, grid_spacing
