# Standard library
from typing import Sequence, Union, Optional, Tuple, Literal
from pathlib import Path

import numpy as np
import jax.numpy as jnp
# Self
from opencadd._typing import PathLike, ArrayLike
from opencadd.const.autodock import AtomType
from opencadd import _exceptions


class GPFFileStructure:

    def __init__(
            self,
            receptor: PathLike,
            receptor_types: Sequence[AtomType] = (
                    AtomType.A, AtomType.C, AtomType.HD, AtomType.N, AtomType.OA, AtomType.SA
            ),
            ligand_types: Sequence[AtomType] = (
                    AtomType.A, AtomType.C, AtomType.HD, AtomType.N, AtomType.NA, AtomType.OA, AtomType.SA
            ),
            gridcenter: Union[Tuple[float, float, float], Literal["auto"]] = "auto",
            npts: Tuple[int, int, int] = (40, 40, 40),
            spacing: float = 0.375,
            smooth: float = 0.5,
            dielectric: float = -0.1465,
            parameter_file: Optional[PathLike] = None,
            output_filename: Optional[str] = None,
            output_path: Optional[PathLike] = None,
    ):
        """
        Parameters
        ----------
        receptor : PathLike
            Filepath to the PDBQT structure file of the macromolecule.
        output_path : PathLike
            Path to a directory to write the output files in.
        receptor_types : sequence of opencadd.const.autodock.AtomType, optional, default: None
            AutoDock atom types present in the receptor. If not specified (i.e. when set to `None`; default),
            atom types will be extracted from the given PDBQT file.
        ligand_types : sequence of opencadd.const.autodock.AtomType, optional, default: (AtomType.C, AtomType.A, AtomType.HD, AtomType.OA)
            Atom types present in the ligand.
        gridcenter : tuple[float, float, float] or "auto", optional, default: "auto"
            Coordinates (x, y, z) of the center of grid map, in angstroms (Å).
            If set to "auto", AutoGrid automatically centers the grid on the center of macromolecule.
        npts : tuple[int, int, int], optional, default: (40, 40, 40)
            Number of grid points to add to the central grid point, along x-, y- and z-axes,
            respectively. Each value must be an even integer number; when added to the central grid
            point, there will be an odd number of points in each dimension. The number of x-, y and
            z-grid points need not be equal.
        spacing : float, optional, default: 0.375
            The grid-point spacing, i.e. distance between two grid points in angstroms (Å).
            Grid points are orthogonal and uniformly spaced in AutoDock, i.e. this value is used for
            all three dimensions.
        smooth : float, optional, default: 0.5
            Smoothing parameter for the pairwise atomic affinity potentials (both van der Waals
            and hydrogen bonds), in angstroms (Å).
            For AutoDock4, the force field has been optimized for a value of 0.5 Å.
        dielectric : float, optional, default: -0.1465
            Dielectric function flag: if negative, AutoGrid will use distance-dependent dielectric
            of Mehler and Solmajer; if the float is positive, AutoGrid will use this value as the
            dielectric constant. AutoDock4 has been calibrated to use a value of –0.1465.
        parameter_file : PathLike, optional, default: None
            User-defined atomic parameter file. If not provided, AutoGrid uses internal parameters.
        """
        self.receptor = receptor
        self.receptor_types = receptor_types
        self.ligand_types = ligand_types
        self.gridcenter = gridcenter
        self.npts = npts
        self.spacing = spacing
        self.smooth = smooth
        self.dielectric = dielectric
        self.parameter_file = parameter_file

        self._output_filename = "temp"
        self.output_path = output_path
        self.output_filename = output_filename
        return

    @property
    def receptor(self) -> Path:
        return self._receptor

    @receptor.setter
    def receptor(self, value: PathLike):
        self._receptor = self._set_filepath(value)
        return

    @property
    def receptor_types(self) -> Tuple[AtomType, ...]:
        return self._receptor_types

    @receptor_types.setter
    def receptor_types(self, value: Union[AtomType, Sequence[AtomType]]):
        self._receptor_types = self._set_atom_types(value)
        return

    @property
    def ligand_types(self):
        return self._ligand_types

    @ligand_types.setter
    def ligand_types(self, value):
        self._ligand_types = self._set_atom_types(value)
        return

    @property
    def gridcenter(self):
        return self._gridcenter

    @gridcenter.setter
    def gridcenter(self, value):
        if isinstance(value, str):
            if value != "auto":
                raise ValueError()
            self._gridcenter = value
            return
        center = np.asarray(value)
        _exceptions.raise_array(
            parent_name=self.__class__.__name__,
            param_name="gridcenter",
            array=center,
            ndim_eq=1,
            size_eq=3,
            dtype=(np.integer, np.floating)
        )
        self._gridcenter = center
        return

    @property
    def npts(self) -> np.ndarray:
        return self._npts

    @npts.setter
    def npts(self, value):
        npts = np.asarray(value)
        _exceptions.raise_array(
            parent_name=self.__class__.__name__,
            param_name="npts",
            array=npts,
            ndim_eq=1,
            size_eq=3,
            dtype=np.integer
        )
        if np.any(npts % 2 != 0):
            raise ValueError()
        self._npts = npts
        return

    @property
    def spacing(self):
        return self._spacing

    @spacing.setter
    def spacing(self, value):
        _exceptions.check_number(value, dtypes="real", gt=0)
        self._spacing = value
        return

    @property
    def smooth(self):
        return self._smooth

    @smooth.setter
    def smooth(self, value):
        _exceptions.check_number(value, dtypes="real", ge=0)
        self._smooth = value
        return

    @property
    def dielectric(self):
        return self._dielectric

    @dielectric.setter
    def dielectric(self, value):
        _exceptions.check_number(value, dtypes="real")
        self._dielectric = value
        return

    @property
    def parameter_file(self):
        return self._parameter_file

    @parameter_file.setter
    def parameter_file(self, value):
        self._parameter_file = None if value is None else self._set_filepath(value)
        return

    @property
    def output_path(self):
        return self._output_path

    @output_path.setter
    def output_path(self, value):
        self._output_path = self.receptor.parent if value is None else self._set_filepath(value)
        self._set_output_paths()
        return

    @property
    def output_filename(self):
        return self._output_filename

    @output_filename.setter
    def output_filename(self, value):
        if value is None:
            self._output_filename = self.receptor.name
        elif isinstance(value, str):
            if " " in value:
                raise ValueError("Path cannot contain spaces.")
            if not value.isascii():
                raise ValueError("Path contains non-ASCII characters.")
            self._output_filename = value
        else:
            raise TypeError()
        self._set_output_paths()
        return

    @property
    def gridfld(self):
        return self._gridfld

    @property
    def xyz(self):
        return self._xyz
    
    @property
    def elecmap(self):
        return self._elecmap
    
    @property
    def dsolvmap(self):
        return self._dsolvmap
    
    @property
    def ligand_maps(self):
        return self._ligand_maps

    def _set_output_paths(self):
        path_common = self.output_path / self.output_filename
        self._gridfld = path_common.with_suffix(".maps.fld")
        self._xyz = path_common.with_suffix(".maps.xyz")
        self._elecmap = path_common.with_suffix(".e.map")
        self._dsolvmap = path_common.with_suffix(".d.map")
        self._ligand_maps = tuple(
            path_common.with_suffix(f'.{ligand_type.name}.map') for ligand_type in self.ligand_types
        )
        return

    @staticmethod
    def _set_atom_types(types: Union[AtomType, Sequence[AtomType]]) -> Tuple[AtomType, ...]:
        if isinstance(types, AtomType):
            return tuple([types])
        if isinstance(types, ArrayLike):
            for atom_type in types:
                if not isinstance(atom_type, AtomType):
                    raise TypeError()
            return tuple(types)
        raise TypeError()

    def _set_filepath(self, filepath: PathLike):
        path = Path(filepath).resolve()
        self._check_filepath(path)
        return path

    @staticmethod
    def _check_filepath(filepath: Path):
        path_str = str(filepath)
        if " " in path_str:
            raise ValueError("Path cannot contain spaces.")
        if not path_str.isascii():
            raise ValueError("Path contains non-ASCII characters.")


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

    if not receptor_filepath.is_file() or receptor_filepath.suffix.lower() != ".pdbqt":
        raise ValueError(f"No PDBQT file found at: {receptor_filepath}")
    if output_path is None:
        output_path = receptor_filepath.parent
    else:
        output_path = Path(output_path).absolute()

        output_path.mkdir(parents=True, exist_ok=True)

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


def write(gpf_file: GPFFileStructure):
    # Generate the file content.
    # It is recommended by AutoDock to generate the gpf file in this exact order.
    file_content: str = ""
    if gpf_file.parameter_file is not None:
        file_content += f"parameter_file {gpf_file.parameter_file}\n"
    file_content += (
        f"npts {' '.join(gpf_file.npts.astype(str))}\n"
        f"gridfld {gpf_file.gridfld}\n"
        f"spacing {gpf_file.spacing}\n"
        f"receptor_types {' '.join(receptor_type.name for receptor_type in gpf_file.receptor_types)}\n"
        f"ligand_types {' '.join(ligand_type.name for ligand_type in gpf_file.ligand_types)}\n"
        f"receptor {gpf_file.receptor}\n"
        f"gridcenter {' '.join(gpf_file.gridcenter.astype(str))}\n"
        f"smooth {gpf_file.smooth}\n"
    )
    for ligand_map in gpf_file.ligand_maps:
        file_content += f"map {ligand_map}\n"
    file_content += (
        f"elecmap {gpf_file.elecmap}\n"
        f"dsolvmap {gpf_file.dsolvmap}\n"
        f"dielectric {gpf_file.dielectric}"
    )
    return file_content

