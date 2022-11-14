"""
T2F-Pharm (Truly Target-Focused Pharmacophore) model.
"""


# Standard library
from typing import Sequence, Union, Optional, Literal, Tuple
from pathlib import Path
# Self
from opencadd.docking import autogrid
from opencadd.io import pdbqt


class T2FPharm:
    """
    Truly Target Focused (T2F) pharmacophore model.
    """

    def __init__(
            self,
            filepath_target_structure: Union[str, Path],
            path_output: Optional[Union[str, Path]] = None,
    ):
        """
        Parameters
        ----------
        filepath_target_structure : str | pathlib.Path
            Path to the input PDBQT file containing the target structure.
        path_output : str | pathlib.Path, Optional, default: None
            An output path to store the results. If not provided, the output files will be
            stored in the same folder as the input file. If a non-existing path is given,
            a new directory will be created with all necessary parent directories.
        """

        # Verify the input file is an existing PDBQT file and assign output path
        filepath_target = Path(filepath_target_structure)
        if not filepath_target.is_file():
            raise ValueError("Input filepath does not point to a file.")
        try:
            df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=filepath_target)
        except Exception as e:
            raise ValueError("Could not parse the input PDBQT file.") from e
        if path_output is None:
            path_output = filepath_target.parent / filepath_target.stem
        else:
            path_output = Path(path_output)
            path_output.mkdir(parents=True, exist_ok=True)

        self._filepath_target = filepath_target
        self._path_output = path_output
        return

    def calculate_energy_grid(
            self,
            probe_types: Sequence[Literal] = ("A", "C", "HD", "OA"),
            coordinates_pocket_center: Union[Tuple[float, float, float], str] = "auto",
            dimensions_pocket: Tuple[float, float, float] = (16, 16, 16),
            spacing_grid: float = 0.6,
            npts: Optional[Tuple[int, int, int]] = None,
    ):
        """

        Parameters
        ----------
        probe_types : Sequence
        coordinates_pocket_center
        dimensions_pocket
        spacing_grid
        npts

        Returns
        -------

        """
        if npts is None:
            npts = autogrid.calculate_npts(
                dimensions_pocket=dimensions_pocket,
                spacing=spacing_grid
            )

        grid = autogrid.routine_run_autogrid(
            receptor=self._filepath_target,
            gridcenter=coordinates_pocket_center,
            npts=npts,
            spacing=spacing_grid,
            ligand_types=probe_types,
            path_output=self._path_output / self._filepath_target.stem
        )
        return grid

