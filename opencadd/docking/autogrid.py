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


import os
from typing import Sequence, Union, Optional, NoReturn
from pathlib import Path


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


def extract_energies_from_map_file(
        filepath: Path
):
    """
    Extract the calculated grid point energies from a '.map' output file into a list.

    Parameters
    ----------
    filepath : pathlib.Path
        Filepath of a map file outputted by AutoGrid.

    Returns
    -------
    list
        A 1-dimensional list with a length equal to the number of grid points. Each element is
        a float, corresponding to the calculated energy value at a grid point. The grid points
        are ordered according to the nested loops $z(y(x))$, so the $x$-coordinate is changing
        fastest. Therefore, the points are ordered according to their coordinates as follows:
        [(0,0,0), (0,0,1), ..., (0,0,nx), (0,1,0), (0,1,1), ..., (0,ny,nx), (1,0,0), (1,0,1), ...,
        (nz,ny,nx)]
    """
    with open(filepath, "r") as f:
        lines = f.readlines()
    # The first 6 lines are headers. The energy values start at line 7,
    # and are written one per line, until the end of file.
    return list(map(float, lines[6:]))
