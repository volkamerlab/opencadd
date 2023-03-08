"""
API for the AutoGrid4 program.

Functions and routines to communicate with the AutoGrid4 program via shell command executions,
and get the results as numpy arrays.

References
----------
https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf
https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf
https://www.csb.yale.edu/userguides/datamanip/autodock/html/Using_AutoDock_305.21.html
"""


# Standard library
from typing import Sequence, Union, Optional, Tuple, Literal
import subprocess
from pathlib import Path
# 3rd-party
import numpy as np
import numpy.typing as npt
# Self
from opencadd._typing import PathLike, ArrayLike
from opencadd.const.autodock import AtomType
from opencadd import spacetime
from opencadd import pocket
import opencadd as oc

_PATH_EXECUTABLE = Path(oc.__file__).parent.resolve()/'_exec'/'autogrid4'


def _from_pdbqt_content(
        content: Sequence[str],
        receptor_types: Sequence[AtomType],
        ligand_types: Sequence[AtomType],
        grid: spacetime.grid.Grid,
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        param_filepath: Optional[Path] = None,
        field_datatype: npt.DTypeLike = np.single,
):
    """
    Run AutoGrid energy calculations and get the results.

    Parameters
    ----------
    content : pathlib.Path
        Path to the PDBQT structure file of the macromolecule.
    ligand_types : Sequence[opencadd.consts.autodock.AtomType]
        Types of ligand atoms, for which interaction energies must be calculated.
        For more information, see `opencadd.consts.autodock.AtomType`.
    grid_center : tuple[float, float, float] | Literal["auto"], Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in the reference frame of the target
        structure, in Ångstrom (Å). If set to "auto", AutoGrid automatically centers the grid
        on the receptor's center of mass.
    grid_npts : Tuple[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    grid_spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two adjacent grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    smooth : float, Optional, default: 0.5
        Smoothing parameter for the pairwise atomic affinity potentials (both van der Waals
        and hydrogen bonds). For AutoDock4, the force field has been optimized for a value of
        0.5 Å.
    dielectric : float, Optional, default: -0.1465
        Dielectric function flag: if negative, AutoGrid will use distance-dependent dielectric
        of Mehler and Solmajer; if the float is positive, AutoGrid will use this value as the
        dielectric constant. AutoDock4 has been calibrated to use a value of –0.1465.
    param_filepath : pathlib.Path, Optional, default: None
        User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.
    output_path: pathlib.Path
        Path to a folder to write the output files in. If not provided, the output files will be
        stored in the same folder as the input file. If a non-existing path is given,
        a new directory will be created with all necessary parent directories.

    Returns
    -------
    field : opencadd.misc.spatial.Field
        Calculated grid-point energies, as a tuple of 1-dimensional arrays containing the
        energy values for each grid point for a specific type of energy. The grid points are
        ordered according to the nested loops z(y(x)), so the x-coordinate is changing fastest.
        The tuple of energy arrays is ordered in the same way as the input `ligand_types`,
        with two additional grids, namely electrostatic potential, and desolvation energy,
        added to the end of the tuple, respectively. The second tuple contains the paths to each
        of the energy map files in the same order.
    """

    num_receptors = len(content)

    center, npts, spacing, slices = _extract_grid_values(grid=grid)

    fields_tensor = np.empty(shape=(num_receptors, *grid.shape, len(ligand_types)+2), dtype=field_datatype)
    temp_path = Path.home() / "opencadd_temp"
    temp_path_pdbqt = temp_path.with_suffix(".pdbqt")
    temp_path_gpf = temp_path.with_suffix(".gpf")
    gpf_file = oc.io.autodock.gpf.GPFFileStructure(
        receptor=temp_path_pdbqt,
        receptor_types=receptor_types,
        ligand_types=ligand_types,
        gridcenter=center,
        npts=npts,
        spacing=spacing,
        smooth=smooth,
        dielectric=dielectric,
        parameter_file=param_filepath,
    )
    with open(temp_path_gpf, "wt") as f:
        f.write(oc.io.autodock.gpf.write(gpf_file))
    for receptor_idx, receptor_filepath in enumerate(content):
        with open(temp_path_pdbqt, "wt") as f:
            f.write(receptor_filepath)
        _submit_job(filepath_gpf=temp_path_gpf)

        paths_maps = gpf_file.ligand_maps + (gpf_file.elecmap, gpf_file.dsolvmap)
        for idx_map, path_map in enumerate(paths_maps):
            with open(path_map) as f:
                lines = f.read().splitlines()
            fields_tensor[receptor_idx, ..., idx_map] = np.array(
                lines[6:], dtype=field_datatype
            ).reshape(npts + 1, order="F")[slices]

    return oc.spacetime.field.from_tensor_grid(
        tensor=fields_tensor,
        grid=grid,
        names=(*(ligand_type.name for ligand_type in ligand_types), "e", "d")
    )


def _extract_grid_values(grid: spacetime.grid.Grid):
    if grid.dimension != 3:
        raise ValueError(f"AutoGrid only accepts 3D grids, but the input grid had {grid.dimension} dimensions.")
    if not np.allclose(grid.spacings, grid.spacings[0]):
        raise ValueError("AutoGrid only accepts grids with equal spacing in all dimensions.")
    is_odd = grid.shape % 2
    if np.all(is_odd):
        return grid.center, grid.shape - 1, grid.spacings[0], slice(None)
    npts = np.where(is_odd, grid.shape - 1, grid.shape)
    center = grid.coordinates[tuple(grid.shape // 2)]
    slices = tuple(np.where(is_odd, slice(None), slice(-1)))
    return center, npts, grid.spacings[0], slices


def _submit_job(
        filepath_gpf: Union[PathLike, Sequence[PathLike]],
        output_path: Optional[PathLike] = None,
) -> subprocess.CompletedProcess:
    """
    Run grid energy calculations with AutoGrid4, using the input grid parameter file (GPF).

    Parameters
    ----------
    filepath_gpf : pathlib.Path
        Path to the input Grid Parameter File (GPF), used as an input specification file in
        AutoGrid.
    output_path : pathlib.Path, Optional, default: None
        Filepath to store the output log of the process produced by AutoGrid. If `None`,
        then the log file will be created in the same directory as `filepath_input_gpf`,
        with the same filename but with .GLG extension.

    Returns
    -------
    subprocess.CompletedProcess
        If the process exits with a zero exit code (meaning it was successful),
        a subprocess.CompletedProcess object is returned, containing the console output and
        other attributes of the process. But more importantly, a log file (.glg) is generated in
        the given output path. Calculated grid-point energies are written to respective MAP
        files, as specified in the input GPF file.

    Raises
    ------
    subprocess.CalledProcessError
        If the process exits with a non-zero exit code, an exception will be raised,
        with attributes `returncode`, `stdout` and `stderr`, which hold the exit code,
        console output and error message of the process.
    """
    if output_path is None:
        output_path = filepath_gpf
    process = subprocess.run(
        args=[
            _PATH_EXECUTABLE,
            "-p",
            Path(filepath_gpf).with_suffix('.gpf'),
            "-l",
            Path(output_path).with_suffix('.glg')
        ],
        capture_output=True,
        check=True,
    )
    return process


def calculate_npts(
        grid_size: Tuple[float, float, float],
        grid_spacing: float
) -> Tuple[int, int, int]:
    """
    Calculate the AutoGrid input argument `npts`.

    Parameters
    ----------
    grid_size : tuple[float, float, float]
        Length of the grid along x-, y-, and z-axis, respectively.
    grid_spacing : float
        The same parameter as in AutoGrid, i.e. the grid-point spacing.

    Returns
    -------
    npts : tuple[int, int, int]
        Can be used directly as input `npts` for AutoGrid functions. The values are the smallest
        valid values (i.e. even integers) that are needed to cover the whole cuboid pocket.
        Therefore, in cases where a dimension is not divisible by the spacing value, or the
        resulting value is an odd number, the value will be rounded up to the next even integer.

    Notes
    -----
    The units of values in `grid_size` and `grid_spacing` don't matter in this function,
    as long as they are both in the same units. Notice that in AutoGrid functions, the `spacing`
    argument must be in Ångstrom.

    See Also
    --------
    For more information on AutoGrid parameters `spacing` and `npts`, see the function
    `routine_run` in this module.
    """
    npts_min = np.ceil(np.array(grid_size) / grid_spacing)
    return tuple(np.where(npts_min % 2 == 0, npts_min, npts_min + 1).astype(int))
