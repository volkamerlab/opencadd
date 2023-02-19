# Standard library
from typing import Sequence, Union, Optional, Tuple, Literal
from pathlib import Path
# Self
from opencadd._typing import PathLike, ArrayLike
from opencadd.const.autodock import AtomType


def write(
        receptor_filepath: PathLike,
        output_path: PathLike = None,
        receptor_types: Sequence[AtomType] = None,
        ligand_types: Sequence[AtomType] = (AtomType.C, AtomType.A, AtomType.HD, AtomType.OA),
        grid_center: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
        grid_npts: Tuple[int, int, int] = (40, 40, 40),
        grid_spacing: float = 0.375,
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        parameter_filepath: Optional[PathLike] = None,
        return_paths: bool = True,
        return_content: bool = False,
) -> Optional[Tuple[Optional[dict], Optional[str]]]:
    """
    Create a Grid Parameter File (GPF), used as an input specification file in AutoGrid.

    Parameters
    ----------
    receptor_filepath : PathLike
        Filepath to the PDBQT structure file of the macromolecule.
    output_path : PathLike
        Path to a folder to write the output files in.
    receptor_types : sequence of opencadd.const.autodock.AtomType, optional, default: None
        AutoDock atom types present in the receptor. If not specified (i.e. when set to `None`; default),
        atom types will be extracted from the given PDBQT file.
    ligand_types : sequence of opencadd.const.autodock.AtomType, optional, default: (AtomType.C, AtomType.A, AtomType.HD, AtomType.OA)
        Atom types present in the ligand.
    grid_center : tuple[float, float, float] or "auto", optional, default: "auto"
        Coordinates (x, y, z) of the center of grid map, in angstroms (Å).
        If set to "auto", AutoGrid automatically centers the grid on the center of macromolecule.
    grid_npts : tuple[int, int, int], optional, default: (40, 40, 40)
        Number of grid points to add to the central grid point, along x-, y- and z-axes,
        respectively. Each value must be an even integer number; when added to the central grid
        point, there will be an odd number of points in each dimension. The number of x-, y and
        z-grid points need not be equal.
    grid_spacing : float, optional, default: 0.375
        The grid-point spacing, i.e. distance between two grid points in angstroms (Å).
        Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
        all three dimensions.
    smooth : float, optional, default: 0.5
        Smoothing parameter for the pairwise atomic affinity potentials (both van der Waals
        and hydrogen bonds). For AutoDock4, the force field has been optimized for a value of
        0.5 Å.
    dielectric : float, optional, default: -0.1465
        Dielectric function flag: if negative, AutoGrid will use distance-dependent dielectric
        of Mehler and Solmajer; if the float is positive, AutoGrid will use this value as the
        dielectric constant. AutoDock4 has been calibrated to use a value of –0.1465.
    parameter_filepath : PathLike, optional, default: None
        User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.


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
    # TODO: apparently AutoGrid cannot handle filepaths in the gpf file that have spaces. Using
    #  quotation marks around the filepath, and escaping with \ did not work. Find a solution.
    if " " in str(receptor_filepath):
        raise ValueError("Path cannot contain spaces.")
    if not receptor_filepath.is_file() or receptor_filepath.suffix.lower() != ".pdbqt":
        raise ValueError(f"No PDBQT file found at: {receptor_filepath}")
    if output_path is None:
        output_path = receptor_filepath.parent
    else:
        output_path = Path(output_path).absolute()
        if " " in str(output_path):
            raise ValueError("Path cannot contain spaces.")
        output_path.mkdir(parents=True, exist_ok=True)
    if receptor_types is None:
        # Parse the AutoDock atom types of receptor from PDBQT file
        with open(receptor_filepath, "r") as f:
            lines = f.readlines()
        # Receptor types are stored in the last column of ATOM records
        atom_types = [line.split()[-1] for line in lines if line.startswith("ATOM")]
        unique_atom_types = set(atom_types)
        receptor_types = tuple(AtomType[atom_type] for atom_type in unique_atom_types)
    # Create filepaths for output files.
    path_common = output_path / receptor_filepath.name
    paths = dict(
        gpf=path_common.with_suffix(".gpf"),
        fld=path_common.with_suffix(".maps.fld"),
        xyz=path_common.with_suffix(".maps.xyz"),
        maps=dict(
            e=path_common.with_suffix(".e.map"),
            d=path_common.with_suffix(".d.map"),
            ligands={
                ligand_type.name: path_common.with_suffix(f'.{ligand_type.name}.map')
                for ligand_type in ligand_types
            }
        )
    )
    # Generate the file content.
    # It is recommended by AutoDock to generate the gpf file in this exact order.
    file_content: str = ""
    if parameter_filepath is not None:
        file_content += f"parameter_file {Path(parameter_filepath).absolute()}\n"
    file_content += (
        f"npts {grid_npts[0]} {grid_npts[1]} {grid_npts[2]}\n"
        f"gridfld {paths['fld']}\n"
        f"spacing {grid_spacing}\n"
        f"receptor_types {' '.join(receptor_type.name for receptor_type in receptor_types)}\n"
        f"ligand_types {' '.join(ligand_type.name for ligand_type in ligand_types)}\n"
        f"receptor {receptor_filepath}\n"
        f"gridcenter {grid_center[0]} {grid_center[1]} {grid_center[2]}\n"
        f"smooth {smooth}\n"
    )
    for _, path_ligand_map in paths['maps']['ligands'].items():
        file_content += f"map {path_ligand_map}\n"
    file_content += (
        f"elecmap {paths['maps']['e']}\n"
        f"dsolvmap {paths['maps']['d']}\n"
        f"dielectric {dielectric}"
    )
    # write the file content to pgf file.
    with open(paths['gpf'], "w") as f:
        f.write(file_content)
    if return_paths:
        if return_content:
            return paths, file_content
        return paths
    elif return_content:
        return file_content
    return
