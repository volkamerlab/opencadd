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
import os
from typing import Sequence, Union, Optional, NoReturn
from pathlib import Path
import itertools
# 3rd-party
import numpy as np


def routine_run_autogrid(
        receptor: Path,
        gridcenter: Union[Sequence[float, float, float], str] = "auto",
        npts: Sequence[int, int, int] = (40, 40, 40),
        spacing: float = 0.375,
        receptor_types: str = "A C HD N NA OA SA Cl",
        ligand_types: str = "A C HD OA",
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        parameter_file: Optional[Path] = None,
        path_output: Optional[Path] = None,
):
    """
    Run AutoGrid with given input structure and specifications, and return the grid values.
    The set of input parameters are identical to that of function `create_gpf`.

    Parameters
    ----------
    receptor : pathlib.Path
        Filepath to the PDBQT structure file of the macromolecule.
    gridcenter : Sequence[float, float, float], Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
        If not provided, the keyword "auto" signals AutoGrid to center the grid on the center of
        the macromolecule.
    npts : Sequence[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    receptor_types : str, Optional, default: "A C HD N NA OA SA Cl"
        Atom types present in the receptor, separated by spaces; e.g. for a typical protein, this
        will be, “A C HD N OA SA”. Atom types are one or two letters, and several specialized
        types are used in the AutoDock4.2 forcefield, including: C (aliphatic carbon), A (aromatic
        carbon), HD (hydrogen that donates hydrogen bond), OA (oxygen that accepts hydrogen
        bond), N (nitrogen that doesn’t accept hydrogen bonds), SA (sulfur that accepts hydrogen
        bonds).
    ligand_types : str, Optional, default: "A C HD OA"
        Atom types present in the ligand, separated by spaces, such as “A C HD N NA OA SA”.
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
    numpy.ndarray, numpy.ndarray
    """
    # 1. Create GPF config file for AutoGrid.
    (
        path_gpf,
        path_gridfld,
        paths_ligand_type_maps,
        path_electrostatic_map,
        path_desolvation_map
    ) = create_gpf(
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
        filepath_output_glg=path_output/receptor.stem,
    )
    # 3. Extract calculated grid-point energies from map files.
    grid = np.empty(
        shape=(npts[2]+1, npts[1]+1, npts[0]+1, len(paths_ligand_type_maps)+2),
        dtype=np.float32
    )
    for idx, filepath_map in enumerate(
            paths_ligand_type_maps + path_electrostatic_map + path_desolvation_map
    ):
        grid[..., idx] = create_grid_tensor_from_map_file(filepath=filepath_map)
    # 4. Calculate coordinates of the grid points in reference frame of the target.
    coordinates = calculate_grid_point_coordinates(
        gridcenter=gridcenter,
        npts=npts,
        spacing=spacing
    )
    return grid, coordinates


def create_gpf(
        receptor: Path,
        gridcenter: Union[Sequence[float, float, float], str] = "auto",
        npts: Sequence[int, int, int] = (40, 40, 40),
        spacing: float = 0.375,
        receptor_types: str = "A C HD N NA OA SA Cl",
        ligand_types: str = "A C HD OA",
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        parameter_file: Optional[Path] = None,
        path_output: Optional[Path] = None,
) -> tuple[Path, Path, list[Path], Path, Path]:
    """
    Create a Grid Parameter File (GPF), which is used as an input specification file in AutoGrid.

    Parameters
    ----------
    receptor : pathlib.Path
        Filepath to the PDBQT structure file of the macromolecule.
    gridcenter : Sequence[float, float, float], Optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in Ångstrom (Å).
        If not provided, the keyword "auto" signals AutoGrid to center the grid on the center of
        the macromolecule.
    npts : Sequence[int, int, int], Optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    spacing : float, Optional, default: 0.375
        The grid-point spacing, i.e. distance between two grid points in Ångstrom (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    receptor_types : str, Optional, default: "A C HD N NA OA SA Cl"
        Atom types present in the receptor, separated by spaces; e.g. for a typical protein, this
        will be, “A C HD N OA SA”. Atom types are one or two letters, and several specialized
        types are used in the AutoDock4.2 forcefield, including: C (aliphatic carbon), A (aromatic
        carbon), HD (hydrogen that donates hydrogen bond), OA (oxygen that accepts hydrogen
        bond), N (nitrogen that doesn’t accept hydrogen bonds), SA (sulfur that accepts hydrogen
        bonds).
    ligand_types : str, Optional, default: "A C HD OA"
        Atom types present in the ligand, separated by spaces, such as “A C HD N NA OA SA”.
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
    tuple[pathlib.Path, pathlib.Path, list[pathlib.Path]]
    Filepath to the generated grid parameter file (.GPF), followed by the filepath to the grid
    field file (.FLD), and a list of filepaths to grid map files (.MAP) for each of the provided
    ligand types in the given order.
    """
    # create filepaths for output files
    path_common = receptor if path_output is None else path_output / receptor.name
    path_gpf, path_gridfld, path_electrostatic_map, path_desolvation_map = (
        path_common.with_suffix(ext) for ext in (".gpf", ".maps.fld", ".e.map", ".d.map")
    )
    paths_ligand_type_maps = [
        path_common.with_suffix(f'.{ligand_type}.map') for ligand_type in ligand_types.split()
    ]
    # Generate the file content.
    # It is recommended by AutoDock to generate the gpf file in this exact order.
    file_content: str = ""
    if parameter_file is not None:
        file_content += f"parameter_file {parameter_file}\n"
    file_content += (
        f"npts {npts[0]} {npts[1]} {npts[2]}\n"
        f"gridfld {path_gridfld}\n"
        f"spacing {spacing}\n"
        f"receptor_types {receptor_types}\n"
        f"ligand_types {ligand_types}\n"
        f"receptor {receptor}\n"
        f"gridcenter {gridcenter}\n"
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
        path_gridfld,
        paths_ligand_type_maps,
        path_electrostatic_map,
        path_desolvation_map
    )


def submit_job(
        filepath_input_gpf: Path,
        filepath_output_glg: Path,
        filepath_autogrid: Path = "/home/michele/Software/autodock4/autogrid4"
) -> NoReturn:
    """
    Initiate grid energy calculations using an input grid parameter file (GPF).

    Parameters
    ----------
    filepath_input_gpf
    filepath_output_glg
    filepath_autogrid

    Returns
    -------
    None
        A log file (.glg) is generated in the given output path.
        Calculated grid-point energies are written to respective '.map' files, as specified in the
        input gpf file.
    """
    # TODO: This function relies on manual installation of AutoGrid4 program and inputting its
    #  filepath here. Find a better option!
    os.system(
        f"{filepath_autogrid} "
        f"-p {filepath_input_gpf.with_suffix('.gpf')} -l {filepath_output_glg.with_suffix('.glg')}"
    )
    return


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
    numpy.ndarray
        A 3-dimensional array of shape (nz, ny, nx), where each point can be indexed by its
        coordinates (normalized to units of grid point spacing). Each element is a float,
        corresponding to the calculated energy value at a grid point. The grid points
        are ordered according to the nested loops $z(y(x))$, so the $x$-coordinate is changing
        fastest. The coordinate system is right-handed.

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
    with open(filepath, "r") as f:
        lines = f.readlines()
    num_points_per_axis = np.array(list(map(float, lines[4].split()[1:]))) + 1
    return np.array(map(float, lines[6:])).reshape(tuple(reversed(num_points_per_axis)))


def calculate_grid_point_coordinates(
        gridcenter: Sequence[float, float, float],
        npts: Sequence[int, int, int],
        spacing: float = 0.375,
):
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
        A 2-dimensional array of shape (n, 3), containing the (x,y,z)-coordinates of n grid
        points.
    """
    origin = np.array(gridcenter) - spacing * np.array(npts) / 2
    num_grid_points_per_axis = np.array(npts) + 1
    coordinates_unit_vectors = itertools.product(
        range(num_grid_points_per_axis[2]),
        range(num_grid_points_per_axis[1]),
        range(num_grid_points_per_axis[0])
    )
    return np.flip(np.array(list(coordinates_unit_vectors)), axis=1) * spacing + origin