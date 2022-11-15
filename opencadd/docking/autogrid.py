"""
API for the AutoGrid4 program.

This module contains functions and routines to communicate with the AutoGrid4 program via
shell command executions.

References
----------
https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf
https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf
https://www.csb.yale.edu/userguides/datamanip/autodock/html/Using_AutoDock_305.21.html
"""


# Standard library
from typing import Sequence, Union, Optional, NoReturn, Tuple, Literal
import subprocess
from pathlib import Path
import itertools
# 3rd-party
import numpy as np
# Self
from opencadd.io.pdbqt import extract_autodock_atom_types_from_pdbqt, TYPE_AUTODOCK_ATOM_TYPE


def routine_run_autogrid(
        receptor: Path,
        gridcenter: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
        npts: Tuple[int, int, int] = (40, 40, 40),
        spacing: float = 0.375,
        ligand_types: Sequence[TYPE_AUTODOCK_ATOM_TYPE] = ("A", "C", "HD", "OA"),
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        parameter_file: Optional[Path] = None,
        path_output: Optional[Path] = None,
) -> np.ndarray:
    """
    Run AutoGrid with given input structure and specifications, and return the grid values.
    The set of input parameters are identical to that of function `create_gpf`.

    Parameters
    ----------
    receptor : pathlib.Path
        Filepath to the PDBQT structure file of the macromolecule.
    gridcenter : tuple[float, float, float] | Literal["auto"], Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
        If not provided, the keyword "auto" signals AutoGrid to center the grid on the center of
        the macromolecule.
    npts : Sequence[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two adjacent grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    ligand_types : Sequence[str], Optional, default: ("A", "C", "HD", "OA")
        Atom types present in the ligand, for example: A, C, HD, N, NA, OA, SA.
    smooth : float, Optional, default: 0.5
        Smoothing parameter for the pairwise atomic affinity potentials (both van der Waals
        and hydrogen bonds). For AutoDock4, the force field has been optimized for a value of
        0.5 Å.
    dielectric : float, Optional, default: -0.1465
        Dielectric function flag: if negative, AutoGrid will use distance-dependent dielectric
        of Mehler and Solmajer; if the float is positive, AutoGrid will use this value as the
        dielectric constant. AutoDock4 has been calibrated to use a value of –0.1465.
    parameter_file : pathlib.Path, Optional, default: None
        User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.
    path_output: pathlib.Path
        Path to a folder to write the output files in.

    Returns
    -------
    grid : numpy.ndarray[dtype=numpy.float64, ndims=4]
        A 4-dimensional array, where the first three dimensions represent the grid points,
        and the last dimension contains the data for each grid point.

        Each grid point contains energy values for each of the input ligand atom-types (probes),
        in the given order. Thus, with 'nt' being the number of input probes, the first 'nt'
        elements of the last dimension of the grid array correspond to their energy values.
        There are five other elements at the end of the last dimension, namely the electrostatic
        potential and desolvation energy, followed by the x-, y- and z-coordinates of the grid
        point in the reference frame of the target structure.

        The shape of the array is thus (nx, ny, nz, nt + 5), where nx, ny, and nz are the number of
        grid points along the x-, y-, and z-axis, respectively. These are each equal to their
        corresponding `npts` value, plus 1 (due to center point). The array is ordered in
        the same way that the grid points are ordered in space, and thus the individual grid
        points can be indexed using their actual coordinates (in unit vectors), assuming the
        origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
        grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.
    """
    # 1. Create GPF config file for AutoGrid.
    receptor_types = tuple(set(extract_autodock_atom_types_from_pdbqt(filepath_pdbqt=receptor)))
    path_gpf, paths_energy_maps, paths_gridfld_xyz = create_gpf(
        receptor=receptor,
        gridcenter=gridcenter,
        npts=npts,
        spacing=spacing,
        receptor_types=receptor_types,
        ligand_types=ligand_types,
        smooth=smooth,
        dielectric=dielectric,
        parameter_file=parameter_file,
        path_output=path_output,
    )
    # 2. Submit job to AutoGrid.
    submit_job(
        filepath_input_gpf=path_gpf,
        filepath_output_glg=path_output,
    )
    # 3. Extract calculated grid-point energies from map files.
    grid = np.empty(
        shape=(*(np.array(npts)+1), len(paths_energy_maps)+3),
        dtype=np.float64
    )
    for idx, filepath_map in enumerate(paths_energy_maps):
        grid[..., idx] = create_grid_tensor_from_map_file(filepath=filepath_map)
    # 4. Calculate coordinates of the grid points in reference frame of the target.
    if gridcenter == "auto":
        gridcenter, _, _ = extract_grid_params_from_mapfile(filepath_map=paths_energy_maps[0])
    grid[..., -3:] = calculate_grid_point_coordinates(
        gridcenter=gridcenter,
        npts=npts,
        spacing=spacing
    )
    return grid


def create_gpf(
        receptor: Path,
        gridcenter: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
        npts: Tuple[int, int, int] = (40, 40, 40),
        spacing: float = 0.375,
        receptor_types: Sequence[TYPE_AUTODOCK_ATOM_TYPE] = ("A", "C", "HD", "N", "NA", "OA", "SA", "Cl"),
        ligand_types: Sequence[TYPE_AUTODOCK_ATOM_TYPE] = ("A", "C", "HD", "OA"),
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        parameter_file: Optional[Path] = None,
        path_output: Optional[Path] = None,
) -> Tuple[Path, Tuple[Path], Path, Path]:
    """
    Create a Grid Parameter File (GPF), used as an input specification file in AutoGrid.

    Parameters
    ----------
    receptor : pathlib.Path
        Filepath to the PDBQT structure file of the macromolecule.
    gridcenter : tuple[float, float, float] | "auto", Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
        If set to "auto", AutoGrid automatically centers the grid on macromolecule's center.
    npts : tuple[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    receptor_types : Sequence[str], Optional, default: ("A", "C", "HD", "N", "NA", "OA", "SA", "Cl")
        Atom types present in the receptor. for a typical protein, this will be: A, C, HD, N, OA,
        SA. Atom types are one or two letters, and several specialized types are used in the
        AutoDock4.2 forcefield, including: C (aliphatic carbon), A (aromatic carbon),
        HD (hydrogen that donates hydrogen bond), OA (oxygen that accepts hydrogen bond),
        N (nitrogen that doesn’t accept hydrogen bonds), SA (sulfur that accepts hydrogen bonds).
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
    parameter_file : pathlib.Path, Optional, default: None
        User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.
    path_output: pathlib.Path
        Path to a folder to write the output files in.

    Returns
    -------
    tuple[pathlib.Path, tuple[pathlib.Path], tuple[pathlib.Path, pathlib.Path]]
    Filepath to the generated grid parameter file (.GPF), followed by a tuple of paths to grid
    map files (.MAP), followed by a tuple of paths to the grid field file (.FLD) and .XYZ file,
    respectively.
    For each input ligand-type there will be a filepath to the corresponding .MAP
    file in the tuple of .MAP files, in the same order as inputted. In addition, the two last
    elements of this tuple correspond to electrostatic and desolvation map files that are always
    generated.
    """
    # Create filepaths for output files.
    path_common = receptor if path_output is None else path_output / receptor.name
    path_gpf, path_gridfld, path_xyz, path_electrostatic_map, path_desolvation_map = (
        path_common.with_suffix(ext) for ext in (".gpf", ".maps.fld", ".maps.xyz", ".e.map", ".d.map")
    )
    paths_ligand_type_maps = [
        path_common.with_suffix(f'.{ligand_type}.map') for ligand_type in ligand_types
    ]
    # Generate the file content.
    # It is recommended by AutoDock to generate the gpf file in this exact order.
    # TODO: apparently AutoGrid cannot handle filepaths in the gpf file that have spaces. Using
    #  quotation marks around the filepath, and escaping with \ did not work. Find a solution.
    file_content: str = ""
    if parameter_file is not None:
        file_content += f"parameter_file {parameter_file}\n"
    file_content += (
        f"npts {npts[0]} {npts[1]} {npts[2]}\n"
        f"gridfld {path_gridfld}\n"
        f"spacing {spacing}\n"
        f"receptor_types {' '.join(receptor_types)}\n"
        f"ligand_types {' '.join(ligand_types)}\n"
        f"receptor {receptor}\n"
        f"gridcenter {gridcenter[0]} {gridcenter[1]} {gridcenter[2]}\n"
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
    return (
        path_gpf,
        tuple(paths_ligand_type_maps + [path_electrostatic_map, path_desolvation_map]),
        (path_gridfld, path_xyz),
    )


def submit_job(
        filepath_input_gpf: Path,
        filepath_output_glg: Optional[Path] = None,
) -> subprocess.CompletedProcess:
    """
    Run grid energy calculations with AutoGrid4, using the input grid parameter file (GPF).

    Parameters
    ----------
    filepath_input_gpf : pathlib.Path
        Path to the input Grid Parameter File (.gpf), used as an input specification file in
        AutoGrid.
    filepath_output_glg : pathlib.Path, Optional, default: None
        Filepath to store the output log of the process produced by AutoGrid. If `None`,
        then the log file will be created in the same directory as `filepath_input_gpf`,
        with the same filename but with '.glg' extension.

    Returns
    -------
    subprocess.CompletedProcess
        If the process exits with a zero exit code (meaning it was successful),
        a subprocess.CompletedProcess object is returned, containing the console output and
        other attributes of the process. But more importan
        A log file (.glg) is generated in the given output path.
        Calculated grid-point energies are written to respective '.map' files, as specified in the
        input gpf file.

    Raises
    ------
    subprocess.CalledProcessError
        If the process exits with a non-zero exit code, an exception will be raised,
        with attributes `returncode`, `stdout` and `stderr`, which hold the exit code,
        console output and error message of the process.
    """
    if filepath_output_glg is None:
        filepath_output_glg = filepath_input_gpf
    process = subprocess.run(
        args=[
            Path(__file__).parent.resolve()/'autogrid4',
            "-p",
            filepath_input_gpf.with_suffix('.gpf'),
            "-l",
            filepath_output_glg.with_suffix('.glg')
        ],
        capture_output=True,
        check=True,
    )
    return process


def create_grid_tensor_from_map_file(
        filepath: Path
):
    """
    Extract the calculated grid point energies from a '.map' output file into a 3-dimensional
    array representing the grid.

    Parameters
    ----------
    filepath : pathlib.Path
        Filepath of a map file outputted by AutoGrid.

    Returns
    -------
    grid : numpy.ndarray[dtype=numpy.float64, ndims=3]
        A 3-dimensional array containing the calculated energy values for each grid point.
        The shape of the array is (nx, ny, nz), where nx, ny, and nz are the number of grid
        points along the x-, y-, and z-axis, respectively. These are each equal to the
        corresponding input `npts` value, plus 1 (due to center point). The array is ordered in
        the same way that the grid points are ordered in space, and thus the individual grid
        points can be indexed using their actual coordinates (in unit vectors), assuming the
        origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
        grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.

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
    The grid points are ordered according to the nested loops z(y(x)), so the x-coordinate
    is changing fastest. The coordinate system is right-handed.
    The header `NELEMENTS` is the same as the input parameter `npts`, defined in function
    `create_gpf`.
    """
    with open(filepath, "r") as f:
        lines = f.readlines()
    num_points_per_axis = np.array(list(map(float, lines[4].split()[1:]))).astype(int) + 1
    return np.array(lines[6:]).astype(np.float64).reshape(tuple(num_points_per_axis), order="F")


def calculate_grid_point_coordinates(
        gridcenter: Tuple[float, float, float],
        npts: Tuple[int, int, int],
        spacing: float = 0.375,
) -> np.ndarray:
    """
    Calculate coordinates of the grid points, in the reference frame of the target structure.
    These are calculated from the same input parameters used in function `create_gpf` to create
    the configuration file for AutoGrid.

    Parameters
    ----------
    gridcenter : Sequence[float, float, float], Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
    npts : Sequence[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.

    Returns
    -------
    numpy.ndarray
        A 4-dimensional array containing the (x,y,z)-coordinates of each grid point
        in Ångstrom (Å), in the target structure's reference frame.

        The shape of the array is (nx, ny, nz, 3), where nx, ny, and nz are the number of grid
        points along the x-, y-, and z-axis, respectively. These are each equal to the
        corresponding input `npts` value, plus 1 (due to center point). The array is ordered in
        the same way that the grid points are ordered in space, and thus the individual grid
        points can be indexed using their actual coordinates (in unit vectors), assuming the
        origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
        grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.
    """
    num_grid_points_per_axis = np.array(npts) + 1
    coordinates_unit_vectors = np.flip(
        np.array(
            list(
                itertools.product(
                    range(num_grid_points_per_axis[2]),
                    range(num_grid_points_per_axis[1]),
                    range(num_grid_points_per_axis[0])
                )
            )
        ),
        axis=1
    ).reshape((*num_grid_points_per_axis, 3), order="F")
    origin = np.array(gridcenter) - spacing * np.array(npts) / 2
    return coordinates_unit_vectors * spacing + origin


def calculate_npts(dimensions_pocket: Tuple[float, float, float], spacing: float) -> np.ndarray:
    """
    Calculate the AutoGrid input argument `npts`, from given pocket dimensions and grid spacing
    values.

    Parameters
    ----------
    dimensions_pocket : tuple[float, float, float]
        Length of the target structure's binding pocket, along x-, y-, and z-axis, respectively.
    spacing : float
        The same parameter as in AutoGrid, i.e. the grid-point spacing.

    Returns
    -------
    npts : numpy.ndarray
        A 1-dimensional array of size 3, which can be used directly as input `npts` for AutoGrid
        functions. The values are the smallest valid values (i.e. even integers) that are needed to
        cover the whole cuboid pocket. Therefore, in cases where a dimension is not divisible by
        the spacing value, or the resulting value is an odd number, the value will be rounded up to
        the next even integer.

    Notes
    -----
    The units of values in `dimensions_pocket` and `spacing` don't matter in this function,
    as long as they are both in the same units. Notice that in AutoGrid functions, the `spacing`
    atgument must be in Ångstrom.

    See Also
    --------
    For more information on AutoGrid parameters `spacing` and `npts`, see the function
    `routine_run_autogrid` in this module.
    """
    npts_min = np.ceil(np.array(dimensions_pocket) / spacing)
    return np.where(npts_min % 2 == 0, npts_min, npts_min + 1).astype(int)


def extract_grid_params_from_mapfile(
        filepath_map: Path
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float], float]:
    """
    Extract `gridcenter`, `npts` and `spacing` parameters from a .MAP file.

    Parameters
    ----------
    filepath_map : pathlib.Path
        Path to one of the calculated .MAP files.

    Returns
    -------
    gridcenter, npts, spacing : tuple[tuple, tuple, float]

    See Also
    --------
    For a sample of '.map' file contents, see function `create_grid_tensor_from_map_file`.
    """
    with open(filepath_map, "r") as f:
        lines = f.readlines()
    spacing = float(lines[3].split()[1])
    npts = tuple(map(float, lines[4].split()[1:4]))
    gridcenter = tuple(map(float, lines[5].split()[1:4]))
    return gridcenter, npts, spacing
