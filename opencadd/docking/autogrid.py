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
from opencadd.typing import PathLike
from opencadd.io import pdbqt
from opencadd.consts import autodock
from opencadd.misc import spatial


def routine_run(
        receptor_filepaths: Sequence[PathLike],
        ligand_types: Sequence[autodock.AtomType],
        grid_center: Union[Tuple[float, float, float]],
        grid_size: Tuple[float, float, float] = (25, 25, 25),
        grid_npts: Optional[Tuple[int, int, int]] = None,
        grid_spacing: float = 0.375,
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        receptor_types: Optional[Sequence[Sequence[autodock.AtomType]]] = None,
        param_filepath: Optional[Path] = None,
        output_path: Optional[Path] = None,
        field_datatype: npt.DTypeLike = np.single,
) -> spatial.Field:
    """
    Run AutoGrid energy calculations and get the results.

    Parameters
    ----------
    receptor_filepaths : pathlib.Path
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
    if isinstance(receptor_filepaths, PathLike):
        receptor_filepaths = [receptor_filepaths]
    elif not isinstance(receptor_filepaths, npt.ArrayLike):
        raise ValueError("`receptor_filepaths` must be path-like or array-like of path-likes.")
    num_receptors = len(receptor_filepaths)
    if receptor_types is None:
        receptor_types = [
            pdbqt.extract_autodock_atom_types_from_pdbqt(filepath_pdbqt=receptor_filepath)
            for receptor_filepath in receptor_filepaths
        ]
    else:
        receptor_types = np.array(receptor_types)
        if receptor_types.ndim == 1:
            receptor_types = np.tile(
                receptor_types, num_receptors
            ).reshape(num_receptors, -1)

    if grid_npts is None:
        grid_npts = calculate_npts(grid_size=grid_size, grid_spacing=grid_spacing)

    grid_shape = np.array(grid_npts) + 1

    fields_tensor = np.empty(
        shape=(num_receptors, *grid_shape, len(ligand_types)+2),
        dtype=field_datatype
    )

    for receptor_idx, receptor_filepath in enumerate(receptor_filepaths):
        # 1. Create GPF config file for AutoGrid.
        path_gpf, paths_fields, path_gridfld, path_xyz = create_gpf(
            receptor_filepath=receptor_filepath,
            output_path=output_path,
            receptor_types=receptor_types[receptor_idx],
            ligand_types=ligand_types,
            grid_center=grid_center,
            grid_npts=grid_npts,
            grid_spacing=grid_spacing,
            smooth=smooth,
            dielectric=dielectric,
            parameter_filepath=param_filepath,
        )
        # 2. Submit job to AutoGrid.
        submit_job(
            gpf_filepath=path_gpf,
            glg_filepath=output_path,
        )
        # 3. Extract calculated grid-point energies from map files.
        for field_idx, path_field in enumerate(paths_fields):
            fields_tensor[receptor_idx, ..., field_idx] = extract_field_values(
                map_filepath=path_field,
                data_type=field_datatype
            ).reshape(shape=tuple(grid_shape), order="F")

    field = spatial.Field(
        field_tensor=fields_tensor,
        grid_origin=np.array(grid_center) - grid_spacing * np.array(grid_npts) / 2,
        grid_point_spacing=grid_spacing,
    )
    return field


def create_gpf(
        receptor_filepath: PathLike,
        output_path: PathLike,
        receptor_types: Sequence[autodock.AtomType],
        ligand_types: Sequence[autodock.AtomType],
        grid_center: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
        grid_npts: Tuple[int, int, int] = (40, 40, 40),
        grid_spacing: float = 0.375,
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        parameter_filepath: Optional[Path] = None,
) -> Tuple[Path, Tuple[Path], Path, Path]:
    """
    Create a Grid Parameter File (GPF), used as an input specification file in AutoGrid.

    Parameters
    ----------
    receptor_filepath : pathlib.Path
        Filepath to the PDBQT structure file of the macromolecule.
    grid_center : tuple[float, float, float] | "auto", Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
        If set to "auto", AutoGrid automatically centers the grid on macromolecule's center.
    grid_npts : tuple[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    grid_spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    ligand_types : Sequence[str], Optional, default: ("A", "C", "HD", "OA")
        Atom types present in the ligand, such as A, C, HD, N, NA, OA, SA.
    smooth : float, Optional, default: 0.5
        Smoothing parameter for the pairwise atomic affinity potentials (both van der Waals
        and hydrogen bonds). For AutoDock4, the force field has been optimized for a value of
        0.5 Å.
    dielectric : float, Optional, default: -0.1465
        Dielectric function flag: if negative, AutoGrid will use distance-dependent dielectric
        of Mehler and Solmajer; if the float is positive, AutoGrid will use this value as the
        dielectric constant. AutoDock4 has been calibrated to use a value of –0.1465.
    parameter_filepath : pathlib.Path, Optional, default: None
        User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.
    output_path: pathlib.Path
        Path to a folder to write the output files in.

    Returns
    -------
    tuple[pathlib.Path, tuple[pathlib.Path], pathlib.Path, pathlib.Path]
    Filepath to the generated grid parameter file (.GPF), followed by a tuple of paths to grid
    map files (.MAP), followed by a tuple of paths to the grid field file (.FLD) and .XYZ file,
    respectively.
    For each input ligand-type there will be a filepath to the corresponding .MAP
    file in the tuple of .MAP files, in the same order as inputted. In addition, the two last
    elements of this tuple correspond to electrostatic and desolvation map files that are always
    generated.
    """
    receptor_filepath = Path(receptor_filepath).absolute()
    output_path = Path(output_path).absolute()
    if not receptor_filepath.is_file() or receptor_filepath.suffix.lower() != ".pdbqt":
        raise ValueError(f"No PDBQT file found at: {receptor_filepath}")
    # TODO: apparently AutoGrid cannot handle filepaths in the gpf file that have spaces. Using
    #  quotation marks around the filepath, and escaping with \ did not work. Find a solution.
    if " " in str(receptor_filepath) + str(output_path):
        raise ValueError("Paths can not contain spaces.")
    # Create filepaths for output files.
    output_path.mkdir(parents=True, exist_ok=True)
    path_common = output_path / receptor_filepath.name
    path_gpf, path_gridfld, path_xyz, path_electrostatic_map, path_desolvation_map = (
        path_common.with_suffix(ext)
        for ext in (".gpf", ".maps.fld", ".maps.xyz", ".e.map", ".d.map")
    )
    paths_ligand_type_maps = [
        path_common.with_suffix(f'.{ligand_type}.map') for ligand_type in ligand_types
    ]
    paths_fields = tuple(paths_ligand_type_maps + [path_electrostatic_map, path_desolvation_map])
    # Generate the file content.
    # It is recommended by AutoDock to generate the gpf file in this exact order.

    file_content: str = ""
    if parameter_filepath is not None:
        file_content += f"parameter_file {parameter_filepath}\n"
    file_content += (
        f"npts {grid_npts[0]} {grid_npts[1]} {grid_npts[2]}\n"
        f"gridfld {path_gridfld}\n"
        f"spacing {grid_spacing}\n"
        f"receptor_types {' '.join(receptor_type.name for receptor_type in receptor_types)}\n"
        f"ligand_types {' '.join(ligand_type.name for ligand_type in ligand_types)}\n"
        f"receptor {receptor_filepath}\n"
        f"gridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]}\n"
        f"smooth {smooth}\n"
    )
    for path_map in paths_ligand_type_maps:
        file_content += f"map {path_map}\n"
    file_content += (
        f"elecmap {path_electrostatic_map}\n"
        f"dsolvmap {path_desolvation_map}\n"
        f"dielectric {dielectric}"
    )
    # write to the file content to pgf file.
    with open(path_gpf, "w") as f:
        f.write(file_content)
    return path_gpf, paths_fields, path_gridfld, path_xyz


def submit_job(
        gpf_filepath: PathLike,
        glg_filepath: Optional[PathLike] = None,
) -> subprocess.CompletedProcess:
    """
    Run grid energy calculations with AutoGrid4, using the input grid parameter file (GPF).

    Parameters
    ----------
    gpf_filepath : pathlib.Path
        Path to the input Grid Parameter File (GPF), used as an input specification file in
        AutoGrid.
    glg_filepath : pathlib.Path, Optional, default: None
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
    if glg_filepath is None:
        glg_filepath = gpf_filepath
    process = subprocess.run(
        args=[
            Path(__file__).parent.resolve()/'autogrid4',
            "-p",
            Path(gpf_filepath).with_suffix('.gpf'),
            "-l",
            Path(glg_filepath).with_suffix('.glg')
        ],
        capture_output=True,
        check=True,
    )
    return process


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
