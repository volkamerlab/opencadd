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
from opencadd.const import autodock
# from opencadd.mif.abc import IntraMolecularInteractionField
#
#
_PATH_EXECUTABLE = Path(__file__).parent.resolve()/'executables'/'autogrid4'


class AutoGridField(IntraMolecularInteractionField):

    def __init__(
            self,
            receptor_filepaths: Sequence[PathLike],
            grid_center: Union[Tuple[float, float, float]],
            grid_size: Tuple[float, float, float] = (25, 25, 25),
            grid_npts: Optional[Tuple[int, int, int]] = None,
            grid_spacing: float = 0.375,
            ligand_type_hydrogen_bond_donor: autodock.AtomType = autodock.AtomType.HD,
            ligand_type_hydrogen_bond_acceptor: autodock.AtomType = autodock.AtomType.OA,
            ligand_type_hydrophobic: autodock.AtomType = autodock.AtomType.C,
            ligand_type_aromatic: autodock.AtomType = autodock.AtomType.A,
            ligand_types_other: Optional[Sequence[autodock.AtomType]] = (),
            smooth: float = 0.5,
            dielectric: float = -0.1465,
            receptor_types: Optional[Sequence[Sequence[autodock.AtomType]]] = None,
            param_filepath: Optional[Path] = None,
            output_path: Optional[Path] = None,
            field_datatype: npt.DTypeLike = np.single,
    ):
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
        elif not isinstance(receptor_filepaths, ArrayLike):
            raise ValueError("`receptor_filepaths` must be path-like or array-like of path-likes.")
        num_receptors = len(receptor_filepaths)

        if receptor_types is not None:
            receptor_types = np.array(receptor_types)
            if receptor_types.ndim == 1:
                receptor_types = np.tile(
                    receptor_types, num_receptors
                ).reshape(num_receptors, -1)
        if grid_npts is None:
            grid_npts = calculate_npts(grid_size=grid_size, grid_spacing=grid_spacing)
        grid_shape = np.array(grid_npts) + 1
        ligand_types=(ligand_type_hydrogen_bond_donor, ligand_type_hydrogen_bond_acceptor,
                      ligand_type_hydrophobic, ligand_type_aromatic, *ligand_types_other)
        fields_tensor = np.empty(
            shape=(num_receptors, *grid_shape, len(ligand_types) + 2),
            dtype=field_datatype
        )
        for receptor_idx, receptor_filepath in enumerate(receptor_filepaths):
            fields_values, paths_fields, path_gpf, path_gridfld, path_xyz = routine_run(
                receptor_filepath=receptor_filepath,
                ligand_types=ligand_types,
                grid_center=grid_center,
                grid_size=grid_size,
                grid_npts=grid_npts,
                grid_spacing=grid_spacing,
                smooth=smooth,
                dielectric=dielectric,
                receptor_types=receptor_types[receptor_idx] if receptor_types is not None else
                None,
                param_filepath=param_filepath,
                output_path=output_path,
                field_datatype=field_datatype,
            )
            for field_idx, field_values in enumerate(fields_values):
                fields_tensor[receptor_idx, ..., field_idx] = field_values.reshape(
                    tuple(grid_shape),
                    order="F"
                )
        super().__init__(
            field_tensor=fields_tensor,
            field_names=[ligand_type.name for ligand_type in ligand_types] + ["Electro.", "Desolv."],
            grid_origin=np.array(grid_center) - grid_spacing * np.array(grid_npts) / 2,
            grid_point_spacing=grid_spacing,
        )
        for lig_idx, ligand_type in enumerate(ligand_types_other):
            setattr(self, ligand_type.name, self._tensor[..., lig_idx+4])
        return

    @property
    def electrostatic(self):
        """
        Electrostatic potential field. The values are in kcal.mol^-1.e^-1.

        Returns
        -------
        numpy.ndarray
            A 4-dimensional array of shape (n_t, n_x, n_y, n_z), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
        """
        return self._tensor[..., -2]

    @property
    def van_der_waals(self):
        """
        Sub-array of `T2FPharm.grid`, containing only the ligand interaction energies.

        Returns
        -------
        numpy.ndarray
            A 5-dimensional array of shape (n_t, n_x, n_y, n_z, n_l), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
            n_l: number of input ligand types.
        """
        return self._tensor[..., 0:4]

    @property
    def h_bond_donor(self):
        return self._tensor[..., 0]

    @property
    def h_bond_acceptor(self):
        return self._tensor[..., 1]

    @property
    def hydrophobic(self):
        return self._tensor[..., 2]

    @property
    def aromatic(self):
        return self._tensor[..., 3]

    @property
    def desolvation(self):
        """
        Desolvation energy field.

        Returns
        -------
        numpy.ndarray
            A 4-dimensional array of shape (n_t, n_x, n_y, n_z), with
            n_t: number of input protein structures.
            n_x, n_y, n_z: number of grid points along x, y, and z directions.
        """
        return self._tensor[..., -1]


def routine_run(
        receptor_filepath: PathLike,
        ligand_types: Sequence[autodock.AtomType],
        grid_center: Union[Tuple[float, float, float]],
        grid_size: Tuple[float, float, float] = (25, 25, 25),
        grid_npts: Optional[Tuple[int, int, int]] = None,
        grid_spacing: float = 0.375,
        smooth: float = 0.5,
        dielectric: float = -0.1465,
        receptor_types: Optional[Sequence[autodock.AtomType]] = None,
        param_filepath: Optional[PathLike] = None,
        output_path: Optional[PathLike] = None,
        field_datatype: npt.DTypeLike = np.single,
) -> Tuple[Tuple[np.ndarray], Tuple[Path], Path, Path, Path]:
    """
    Run AutoGrid energy calculations and get the results.

    Parameters
    ----------
    receptor_filepath : pathlib.Path
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

    if grid_npts is None:
        grid_npts = calculate_npts(grid_size=grid_size, grid_spacing=grid_spacing)
    # 1. Create GPF config file for AutoGrid.
    path_gpf, paths_fields, path_gridfld, path_xyz = create_gpf(
        receptor_filepath=receptor_filepath,
        output_path=output_path,
        receptor_types=receptor_types,
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
    fields_values = []
    for field_idx, path_field in enumerate(paths_fields):
        fields_values.append(
            extract_field_values(map_filepath=path_field, data_type=field_datatype)
        )
    return tuple(fields_values), paths_fields, path_gpf, path_gridfld, path_xyz


def from_filepath_gpf(
        filepath: Union[PathLike, Sequence[PathLike]],
        output_path: Optional[PathLike] = None,
) -> subprocess.CompletedProcess:
    """
    Run grid energy calculations with AutoGrid4, using the input grid parameter file (GPF).

    Parameters
    ----------
    filepath : pathlib.Path
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
        output_path = filepath
    process = subprocess.run(
        args=[
            _PATH_EXECUTABLE,
            "-p",
            Path(filepath).with_suffix('.gpf'),
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
