
# Standard library
from typing import Sequence, Union, Optional, Tuple, Literal
import subprocess
from pathlib import Path
# 3rd-party
import numpy as np
# Self
from opencadd._typing import PathLike, ArrayLike
from opencadd import spacetime


class AutoDockMAPParsingError(Exception):
    pass


def from_filepath(
        filepath: Union[PathLike, Sequence[PathLike], Sequence[Sequence[PathLike]]],
        data_type: np.dtype = np.single,
        name: Optional[Union[str, Sequence[str]]] = None,
):
    """
    Extract the AutoGrid parameters `gridcenter`, `npts` and `spacing` from a MAP file.
    Extract the calculated grid point energies from a MAP output file.


    Parameters
    ----------
    filepath
    data_type : numpy.dtype, optional, default: numpy.single
        Numpy datatype of the output array. Default is 32-bit float (numpy.single).
    name

    map_filepath : pathlib.Path
        Filepath of a MAP file outputted by AutoGrid. Path to one of the calculated MAP files.

    Returns
    -------

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
    A 1-dimensional array containing the calculated energy values for each grid point.
    The grid points are ordered according to the nested loops z(y(x)), so the x-coordinate
    is changing fastest. The coordinate system is right-handed.

    The header `NELEMENTS` is the same as the input parameter `npts`, defined in function
    `create_gpf`.
    """
    if isinstance(filepath, PathLike):
        filepath = [[filepath]]
    else:
        if not isinstance(filepath, ArrayLike):
            raise TypeError("")
        if np.any([isinstance(filepath_, PathLike) for filepath_ in filepath]):
            if not np.all([isinstance(filepath_, PathLike) for filepath_ in filepath]):
                raise TypeError("Heterogeneous sequence.")
            filepath = [filepath]
        else:
            if not np.all([isinstance(filepath_, ArrayLike) for filepath_ in filepath]):
                raise TypeError("Heterogeneous sequence.")
    lengths = [len(filepath_) for filepath_ in filepath]
    if np.unique(lengths).size != 1:
        raise TypeError("Ragged nested sequences.")

    with open(filepath[0][0], "r") as f:
        lines = f.read().splitlines()

    for line in lines:
        if line.startswith("NELEMENTS"):
            npts = np.array(line.split()[1:4], dtype=np.short)
            break
    else:
        raise AutoDockMAPParsingError()
    grid_shape = npts + 1
    fields = np.empty(shape=(len(filepath), *grid_shape, lengths[0]), dtype=data_type)
    spacings = np.empty(shape=(len(filepath), lengths[0]))
    npts = np.empty(shape=(len(filepath), lengths[0], 3), dtype=np.short)
    centers = np.empty(shape=(len(filepath), lengths[0], 3))
    paths_gpf = np.empty(shape=(len(filepath), lengths[0]), dtype=object)
    paths_fld = np.empty(shape=(len(filepath), lengths[0]), dtype=object)
    paths_pdbqt = np.empty(shape=(len(filepath), lengths[0]), dtype=object)
    for idx_instance, instance in enumerate(filepath):
        for idx_map, map_path in enumerate(instance):
            grid_point_values, tokens = _parse_single_file(filepath=map_path, data_type=data_type)
            fields[idx_instance, ..., idx_map] = grid_point_values.reshape(grid_shape, order="F")
            spacings[idx_instance, idx_map] = tokens["SPACING"]
            npts[idx_instance, idx_map] = tokens["NELEMENTS"]
            centers[idx_instance, idx_map] = tokens["CENTER"]
            paths_gpf[idx_instance, idx_map] = tokens.get("GRID_PARAMETER_FILE")
            paths_fld[idx_instance, idx_map] = tokens.get("GRID_DATA_FILE")
            paths_pdbqt[idx_instance, idx_map] = tokens.get("MACROMOLECULE")
    if np.any(spacings != spacings[0]):
        raise AutoDockMAPParsingError
    if np.any(npts != npts[0]):
        raise AutoDockMAPParsingError
    if np.any(centers != centers[0]):
        raise AutoDockMAPParsingError
    grid = spacetime.grid.from_center_spacing_shape(
        center=centers[0, 0], spacings=spacings[0, 0], shape=grid_shape
    )
    return spacetime.field.from_tensor_grid(tensor=fields, grid=grid)


def _parse_single_file(filepath: PathLike, data_type: np.dtype):
    with open(filepath, "r") as f:
        lines = f.read().splitlines()
    tokens = dict()
    for line_idx, line in enumerate(lines):
        try:
            float(line)
        except ValueError:
            split = line.split()
            if len(split) in [0, 1]:
                _raise_or_warn()
            elif split[0] in ("GRID_PARAMETER_FILE", "GRID_DATA_FILE", "MACROMOLECULE"):
                tokens[split[0]] = Path(split[1])
            elif split[0] == "SPACING":
                tokens[split[0]] = float(split[1])
            elif split[0] == "NELEMENTS":
                tokens[split[0]] = np.array(split[1:4]).astype(np.short)
            elif split[0] == "CENTER":
                tokens[split[0]] = np.array(split[1:4]).astype(data_type)
            else:
                _raise_or_warn()
        else:
            grid_point_values = np.array(lines[line_idx:], dtype=data_type)
            break
    else:
        raise AutoDockMAPParsingError(
            f"No grid point values found in MAP file at {filepath}."
        )
    keys = tokens.keys()
    for token in ("SPACING", "NELEMENTS", "CENTER"):
        if token not in keys:
            raise AutoDockMAPParsingError()
    return grid_point_values, tokens


def _raise_or_warn(msg, raise_level, strictness):
    pass