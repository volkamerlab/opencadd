"""
T2F (Truly Target-Focused) Pharmacophore modeler.

This module contains the `T2Fpharm` class, used for pharmacophore modeling
from protein apo structures.
"""


# Standard library
from typing import Literal, Optional, Sequence, Tuple, List, Union, Any, Callable
import operator
import asyncio
from time import time
# 3rd-party
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from ipywidgets import interact
# Self
from opencadd.chemical import protein
from opencadd import interaction
from opencadd.consts import autodock
from opencadd.spacetime import field, vectorized
from opencadd.visualization import nglview_api
import ipywidgets as widgets
from IPython.display import display


class T2FPharm:
    """
    Truly Target Focused (T2F) pharmacophore modeler.

    Design a pharmacophore from one or several protein apo structures.

    Instantiate the modeler by entering one or several protein structures, defining their
    binding pocket's coordinates and size, and indicating a set of pharmacophore features
    to consider. The modeler automatically calculates interaction energy fields inside the
    binding pocket of each protein structure, for each pharmacophore feature, plus desolvation
    energy and electrostatic potential fields. All data are stored in a single higher-dimensional
    array, allowing for easy and fast data manipulation and visualization. With the methods
    available to the class, a variety of information on the binding pocket can be extracted and
    used to develop a pharmacophore model.



    calculate the interaction energies
    between each protein structure and each atomic probe (representing a pharmacophore feature),
    at regular grid points spanned over the defined binding pocket.

    A 3-dimensional array containing the calculated energy values for each grid point.
        The shape of the array is (nx, ny, nz), where nx, ny, and nz are the number of grid
        points along the x-, y-, and z-axis, respectively. These are each equal to the
        corresponding input `npts` value, plus 1 (due to center point). The array is ordered in
        the same way that the grid points are ordered in space, and thus the individual grid
        points can be indexed using their actual coordinates (in unit vectors), assuming the
        origin is at the edge of the grid with smallest x, y, and z values. That is, indexing
        grid[i, j, k] gives the point located at position (i, j, k) on the actual grid.


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

    Get interaction energies between the target structure and each of the input probe types
        (plus electrostatic and desolvation energies), calculated using AutoGrid4 at every point
        on a regular (i.e. evenly spaced) cuboid grid, placed on (usually the binding pocket of)
        the target structure, as defined by the input arguments.
    """

    def __init__(
            self,
            receptor: protein.Protein,
            fields: interaction.abc.IntraMolecularInteractionField,
            distance_limit: float = 5,
    ):
        """
        Parameters
        ----------
        receptor : opencadd.protein.Protein
            Receptor of the target-focused pharmacophore model.
        fields : opencadd.misc.field.ToxelField
            Toxel fields calculated for the receptor.
        """
        # Declare attributes
        self._receptor: protein.Protein
        self._fields: field.ToxelField
        self._dist_limit_curr: float
        self._dist_atoms_within_limit: np.ndarray
        self._dist_nearest_atoms: np.ndarray
        self._idx_nearest_atoms: np.ndarray
        # Assign attributes
        self._receptor = receptor
        self._fields = fields
        self._dist_limit_curr = distance_limit
        # Distance-related attributes
        shape_distance_data = (
            self._receptor.trajectory_length,
            *self._fields.grid.shape,
            self._receptor.count_atoms
        )
        self._dist_atoms_within_limit = self._receptor.distance_to_atoms_sparse(
            kdtree=self._fields.grid.kdtree,
            max_distance=self._dist_limit_curr
        ).reshape(shape_distance_data)
        nearest_atoms_distances, nearest_atoms_indices = self._receptor.nearest_atoms(
            coords=self._fields.grid.coordinates.reshape(-1, 3)
        )
        self._dist_nearest_atoms = nearest_atoms_distances.reshape(shape_distance_data[:-1])
        self._idx_nearest_atoms = nearest_atoms_indices.reshape(shape_distance_data[:-1])

        self._nglwidget = self._receptor.create_new_ngl_widget()

        self._vacancy: np.ndarray = None
        self._psp_distances: np.ndarray = None
        self._buriedness: np.ndarray = None
        self._distances_to_protein_atoms: np.ndarray = None
        self._hbonding_atoms_count: np.ndarray = None
        return

    @property
    def receptor(self) -> protein.Protein:
        """
        Receptor of the model.
        """
        return self._receptor

    @property
    def fields(self):
        """
        Fields of the model.0
        """
        return self._fields

    @property
    def distances_to_nearest_protein_atom(self) -> np.ndarray:
        return self._dist_nearest_atoms





    @property
    def ngl_widget(self):
        return self._nglwidget

    @property
    def distances_to_protein_atoms(self):
        return (
            self._distances_to_protein_atoms if self._distances_to_protein_atoms is not None
            else self.calculate_distances_to_protein_atoms()
        )

    @property
    def vacancy(self):
        return self._vacancy if self._vacancy is not None else self.calculate_vacancy()

    @property
    def psp_distances(self):
        return self._psp_distances if self._psp_distances is not None else self.calculate_psp_distances()

    @property
    def buriedness(self):
        return self._buriedness if self._buriedness is not None else self.calculate_buriedness()

    @property
    def hbonding_atoms_count(self):
        return self._hbonding_atoms_count if self._hbonding_atoms_count is not None else self.count_hbonding_atoms_in_radius()

    def filter_by_dist_nearest_atom(
            self,
            dist_range: Tuple[float, float],
            invert: bool = False
    ) -> np.ndarray:
        """
        Filter points by their distances to their nearest atoms.

        Parameters
        ----------
        dist_range : tuple[float, float]
            Lower and upper bounds of distance to the nearest atom.
        invert : bool, optional, default: False
            If False (default), the points whose nearest atoms lie within the distance range are
            returned. If True, the points whose nearest atoms lie outside the range are returned.

        Returns
        -------
        boolean_mask : ndarray
        """
        op_lower, op_upper = (operator.le, operator.ge) if invert else (operator.ge, operator.le)
        op_combine = np.logical_or if invert else np.logical_and
        mask = op_combine(
            op_lower(self.distances_to_nearest_protein_atom, dist_range[0]),
            op_upper(self.distances_to_nearest_protein_atom, dist_range[1]),
        )
        return mask

    def __call__(
            self,
            vacancy_max_energy: float = +0.6,
            buriedness_psp_num_dirs: Literal[3, 7, 13] = 7,
            buriedness_psp_max_len: float = 10.0,
            buriedness_psp_min_count: int = 4,
            hbond_max_len: float = 3.0,
            probes_max_energy: Tuple[float, float, float, float] = (-0.6, -0.35, -0.4, -0.4),
            electrostat_pot_exclusion_range: Tuple[float, float] = (-1.0, 1.0),
            min_neighbor_dist_clustering: float = 1.21,
            min_common_neighbor_count_clustering: int = 6,
            min_points_per_cluster_count: int = 15,
    ):
        hba = autodock.AtomType.OA  # hydrogen-bond acceptor probe
        hbd = autodock.AtomType.HD  # hydrogen-bond donor probe
        aliph = autodock.AtomType.C  # aliphatic hydrophobic probe
        arom = autodock.AtomType.A  # aromatic hydrophobic probe
        probes = (hba, hbd, aliph, arom)

        # Vacancy of each grid point as a boolean grid
        vacant = self.vacancy(energy_cutoff=vacancy_max_energy)
        # Buriedness of each vacant point as a boolean array
        buried = self.calculate_buriedness(
            vacancy=vacant,
            psp_distances=self.grid_psp_dist(
                grid_vacancy=vacant,
                num_directions=buriedness_psp_num_dirs,
                max_radius=buriedness_psp_max_len
            ),
            psp_max_length=buriedness_psp_max_len,
            psp_min_count=buriedness_psp_min_count
        )
        # A boolean grid indicating whether a specific probe at a specific grid point has high affinity
        is_high_affinity = self.grid_probes < np.array(probes_max_energy)
        # A boolean grid indicating points that are both high affinity and buried
        high_affinity_n_buried = np.logical_and(
            is_high_affinity,
            buried[..., np.newaxis],
        )
        # Count of H-bond acceptor and donor atoms (from protein) in proximity of each grid point
        count_hba, count_hbd = self.count_hbonding_atoms_in_radius(
            max_len_hbond=hbond_max_len,
        )
        # Boolean grids describing the H-bond environment for each grid point
        has_hba_env = count_hba > 0
        has_hbd_env = count_hbd > 0
        has_hydrophilic_env = np.logical_or(has_hba_env, has_hbd_env)
        has_hydrophobic_env = np.logical_not(has_hydrophilic_env)

        is_hba, is_hbd, is_aliph, is_arom = (
            np.logical_and(high_affinity_n_buried[..., probes.index(probe)], condition)
            for probe, condition in
        zip(probes, [has_hbd_env, has_hba_env, has_hydrophobic_env, has_hydrophobic_env])
        )

        is_buried_in_hydrophilic_env = np.logical_and(buried, has_hydrophilic_env)
        is_pi = np.logical_and(
            is_buried_in_hydrophilic_env,
            self.grid_electrostatic < electrostat_pot_exclusion_range[0],
        )
        is_ni = np.logical_and(
            is_buried_in_hydrophilic_env,
            self.grid_electrostatic > electrostat_pot_exclusion_range[1],
        )
        return

    def calculate_vacancy(
            self,
            source: Union[str, Sequence[str]] = "mindist",
            vacancy_range: Tuple[float, float] = (2, ),
            energy_cutoff: float = +0.6,
            mode: Optional[Literal["max", "min", "avg", "sum"]] = "min",
    ) -> np.ndarray:
        """
        Calculate whether each grid point is vacant, or occupied by a target atom.

        Parameters
        ----------
        energy_cutoff : float, Optional, default: +0.6
            Cutoff value for energy; grid points with energies lower than cutoff are considered
            vacant.
        mode: Literal["max", "min", "avg", "sum"], Optional, default: "min"
            If the energy of more than one ligand type is to be compared, this parameter defines
            how those different energy values must be processed, before comparing with the cutoff.
        ligand_types : Sequence[opencadd.consts.autodock.AtomType], Optional, default: None
            A subset of ligand types that were used to initialize the object, whose energy values
            are to be taken as reference for calculating the vacancy of each grid point. If not
            set to None, then all ligand interaction energies are considered.

        Returns
        -------
        vacancy : numpy.ndarray[dtype=numpy.bool_, shape=T2FPharm.grid.shape[:-1]]
            A 4-dimensional boolean array matching the first four dimensions of `T2FPharm.grid`,
            indicating whether each grid point is vacant (True), or occupied (False).
            Vacant grid points can easily be indexed by `T2FPharm.grid[vacancy]`.
        """
        # The reducing operations corresponding to each `mode`:
        red_fun = {"max": np.max, "min": np.min, "avg": np.mean, "sum": np.sum}
        # Get index of input ligand types
        # if ligand_types is None:
        #     ind = slice(None)
        # else:
        #     ind = np.argwhere(np.expand_dims(ligand_types, axis=-1) == self._probe_types)[:, 1]
        #     # Verify that all input ligand types are valid
        #     if len(ind) != len(ligand_types):
        #         raise ValueError(f"Some of input energies were not calculated.")
        # Reduce the given references using the given operation.
        energy_vals = red_fun[mode](self._fields.van_der_waals, axis=-1)
        # Apply cutoff and return
        self._vacancy = energy_vals < energy_cutoff
        return self._vacancy


    def calculate_buriedness(
            self,
            vacancy: Optional[np.ndarray] = None,
            psp_distances: Optional[np.ndarray] = None,
            psp_max_length: float = 10.0,
            psp_min_count: int = 4,
    ) -> np.ndarray:
        """
        Calculate whether each grid point is buried inside the target structure or not, based on
        counting the number of protein-solvent-protein (PSP) events for each point, and applying
        a cutoff.

        Notice that before using this method, a vacancy grid must be calculated using the method
        `calculate_vacancy_from_energy`.

        Parameters
        ----------
        vacancy : numpy.ndarray

        psp_distances : numpy.ndarray
        psp_max_length : float, Optional, default: 10.0
            Maximum acceptable distance for a PSP event, in Ångstrom (Å).
        psp_min_count : int, Optional, default: 4
            Minimum required number of PSP events for a grid point, in order to count as buried.

        Returns
        -------
        buriedness : numpy.ndarray[dtype=numpy.bool_, shape=grid_vacancy.shape]
            A 4-dimensional array matching the dimensions of the grid,
            indicating whether each grid point is buried (True), or exposed (False).
            Thus, the buried grid points can easily be indexed using boolean indexing with this
            array, i.e.: `grid[buriedness]`.
        """
        if vacancy is None:
            vacancy = self.vacancy
        if psp_distances is None:
            psp_distances = self.psp_distances
        # Count PSP events that are shorter than the given cutoff.
        grid_psp_counts = np.count_nonzero(psp_distances <= psp_max_length, axis=-1)
        buriedness = grid_psp_counts >= psp_min_count
        buriedness[np.logical_not(vacancy)] = False
        return buriedness

    def calculate_distances_to_protein_atoms(self):
        self._distances_to_protein_atoms = self._fields.grid.distance(
            coordinates=self._receptor.trajectory
        )
        return self._distances_to_protein_atoms

    def count_hbonding_atoms_in_radius(
            self,
            max_len_hbond: float = 3.0,
    ) -> np.ndarray:
        proximates = self.distances_to_protein_atoms < max_len_hbond
        self._hbonding_atoms_count = np.empty(
            shape=(*self.distances_to_protein_atoms.shape[:-1], 2),
            dtype=np.byte
        )
        for idx, hbond_role in enumerate(["hbond_acc", "hbond_don"]):
            self._hbonding_atoms_count[..., idx] = np.count_nonzero(
                np.logical_and(self._receptor.atom_data[hbond_role], proximates),
                axis=-1
            )
        return self._hbonding_atoms_count

    @property
    def interactive_view(self):
        view = self._receptor.create_new_ngl_widget()

        def visualizer(
                vacancy_energy_cutoff: float,
                vacancy_calc_mode: [Literal["max", "min", "avg", "sum"]],
                range_psp_len: Tuple[float, float],
                psp_len_min: float,
                psp_len_max: float,
                psp_count_min: int,
                psp_count_max: int,
                energy_type,
                energy_min: float,
                energy_max: float,

        ):
            nglview_api.remove_component_by_name(view=view, name="grid")
            self.calculate_vacancy(energy_cutoff=vacancy_energy_cutoff, mode=vacancy_calc_mode)
            self.calculate_psp_distances()
            mask1 = np.logical_and(
                self.psp_distances >= psp_len_min,
                self.psp_distances <= psp_len_max
            )
            psp_counts = np.count_nonzero(mask1, axis=-1)
            mask2 = np.logical_and(psp_counts >= psp_count_min, psp_counts <= psp_count_max)
            mask3 = np.logical_and(energy_type >= energy_min, energy_type <= energy_max)
            mask = np.logical_and(np.logical_and(mask1, mask2), mask3)[0]
            self.visualize(grid_mask=mask)

        return interact(
            visualizer,
            vacancy_energy_cutoff=(
                np.min(self.field.van_der_waals),
                np.max(self.field.van_der_waals),
                0.1,
            ),
            vacancy_calc_mode=["max", "min", "avg", "sum"],
            psp_len_min=(0, 50, 1),
            psp_len_max=(0, 50, 1),
            psp_count_min=(0, 13, 1),
            psp_count_max=(0, 13, 1),
            energy_type = [("HD", self.field.h_bond_donor), ("HA", self.field.h_bond_acceptor)],
            energy_min=(-2,10,1),
            energy_max=(-2,10,1),
        )


    def visualize(
            self,
            grid_mask = None,
            weights1 = None,
            weights2 = None,
            color_map: str = "bwr",
            opacity: float = 0.8,
    ):
        if grid_mask is None:
            grid_mask = np.ones(shape=self._fields.grid.shape, dtype=np.bool_)
        normalizer = mpl.colors.Normalize()
        mapper = plt.cm.ScalarMappable(
            norm=normalizer,
            cmap=mpl.colormaps[color_map]
        )
        if weights1 is None:
            weights1 = np.repeat(0.5, np.count_nonzero(grid_mask) * 3)
        elif isinstance(weights1, Sequence):
            weights1 = np.tile(weights1, np.count_nonzero(grid_mask) * 3)
        else:
            weights1 = (mapper.to_rgba(weights1.flatten())[..., :3])
        if weights2 is None:
            weights2 = (
                    np.ones(np.count_nonzero(grid_mask)) * self._fields.grid.spacing / 4
            )
        elif isinstance(weights2, (int, float)):
            weights2 = np.ones(np.count_nonzero(grid_mask)) + weights2
        else:
            weights2 = normalizer(weights2).flatten() * (self._fields.grid.spacing - 0.05) + 0.05
        nglview_api.add_spheres(
            view=self._nglwidget,
            coords=self._fields.grid.coordinates[grid_mask],
            colors=weights1,
            radii=weights2,
            opacity=opacity
        )
        return


class Timer:
    def __init__(self, timeout, callback):
        self._timeout = timeout
        self._callback = callback

    async def _job(self):
        await asyncio.sleep(self._timeout)
        self._callback()

    def start(self):
        self._task = asyncio.ensure_future(self._job())

    def cancel(self):
        self._task.cancel()

def throttle(wait):
    """ Decorator that prevents a function from being called
        more than once every wait period. """
    def decorator(fn):
        time_of_last_call = 0
        scheduled, timer = False, None
        new_args, new_kwargs = None, None
        def throttled(*args, **kwargs):
            nonlocal new_args, new_kwargs, time_of_last_call, scheduled, timer
            def call_it():
                nonlocal new_args, new_kwargs, time_of_last_call, scheduled, timer
                time_of_last_call = time()
                fn(*new_args, **new_kwargs)
                scheduled = False
            time_since_last_call = time() - time_of_last_call
            new_args, new_kwargs = args, kwargs
            if not scheduled:
                scheduled = True
                new_wait = max(0, wait - time_since_last_call)
                timer = Timer(new_wait, call_it)
                timer.start()
        return throttled
    return decorator


def debounce(wait):
    """ Decorator that will postpone a function's
        execution until after `wait` seconds
        have elapsed since the last time it was invoked. """
    def decorator(fn):
        timer = None
        def debounced(*args, **kwargs):
            nonlocal timer
            def call_it():
                fn(*args, **kwargs)
            if timer is not None:
                timer.cancel()
            timer = Timer(wait, call_it)
            timer.start()
        return debounced
    return decorator


class T2FPharmWidget:

    def __init__(
            self,
            modeler: T2FPharm
    ):
        # DECLARE ATTRIBUTES
        self._t2fpharm: T2FPharm
        # Widgets for: BINDING SITE REFINEMENT -> VACANCY
        self._widget_refine_vacancy_source_buttons: List[widgets.ToggleButton]
        self._widget_refine_vacancy_mode_buttons: widgets.ToggleButtons
        self._widget_refine_vacancy_filter_buttons: widgets.ToggleButtons
        self._widget_refine_vacancy_range_slider: widgets.FloatRangeSlider
        self._widget_refine_vacancy_range_unit:  widgets.Label
        # Widgets for: BINDING SITE REFINEMENT -> BURIEDNESS
        self._widget_refine_buriedness_source_buttons: List[widgets.ToggleButton]
        self._widget_refine_buriedness_mode_buttons: widgets.ToggleButtons
        self._widget_refine_buriedness_filter_buttons: widgets.ToggleButtons
        self._widget_refine_buriedness_range_slider: widgets.FloatRangeSlider
        self._widget_refine_buriedness_range_unit:  widgets.Label

        # Assign attributes
        self._t2fpharm = modeler

        self._vacancy_data = np.concatenate(
            (
                self._t2fpharm.fields.tensor,
                self._t2fpharm.distances_to_nearest_protein_atom[..., np.newaxis]
            ),
            axis=-1
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY -> SOURCE -> TOGGLE BUTTONS
        # These are toggle buttons for selecting the sources of vacancy calculation.
        # First add the buttons that are always present, regardless of the fields,
        # then add buttons for each field:
        self._widget_refine_vacancy_source_buttons = self._create_toggle_buttons_source(
            labels=["Distance"]+list(self._t2fpharm.fields.field_names),
            tooltips=["Distance of grid points to protein atoms."]+[
                f"Field {field_name}" for field_name in self._t2fpharm.fields.field_names
            ],
            observer=self._on_value_change_vacancy_source
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY -> MODE -> TOGGLE BUTTONS
        # These are toggle buttons for selecting the mode of vacancy calculation.
        self._widget_refine_vacancy_mode_buttons = self._create_toggle_buttons_mode(
            observer=self._on_value_change_vacancy_mode
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY -> FILTER -> TOGGLE BUTTONS
        # These are toggle buttons for deciding whether the given value range should be included
        # or excluded.
        self._widget_refine_vacancy_filter_buttons = self._create_toggle_buttons_filter(
            observer=self._on_value_change_vacancy_filter
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY -> RANGE -> SLIDER
        # This is a range slider for selecting the range of values for vacancy calculation:
        self._widget_refine_vacancy_range_slider = self._create_range_slider(
            observer=self._on_value_change_vacancy_range
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY -> RANGE -> UNIT
        # This is a label that shows the unit of the selected vacancy source values:
        self._widget_refine_vacancy_range_unit = widgets.Label(
            value="",
            layout=widgets.Layout(margin="0px 0px 0px 0px", align="left")
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY
        # This is a container Box for all vacancy widgets above, plus added labels, making the
        # vacancy control panel:
        vacancy_control_panel = self._create_control_panel(
            controllers=[
                widgets.HBox(self._widget_refine_vacancy_source_buttons),
                widgets.Box(
                    children=[self._widget_refine_vacancy_mode_buttons],
                    layout=widgets.Layout(
                        width='330px',
                        height='',
                        flex_flow='row',
                        display='flex'
                    )
                ),
                self._widget_refine_vacancy_filter_buttons,
                widgets.HBox(
                    [
                        self._widget_refine_vacancy_range_slider,
                        self._widget_refine_vacancy_range_unit
                    ]
                ),
            ],
            header_name="Vacancy"
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> BURIEDNESS -> SOURCE -> TOGGLE BUTTON
        # These are toggle buttons for selecting the sources of buriedness calculation.
        self._widget_refine_buriedness_source_buttons = self._create_toggle_buttons_source(
            labels=["PSP Lengths"],
            tooltips=["Length of PSP events in all 13 cubic directions."],
            observer=self._on_value_change_buriedness_source
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> BURIEDNESS -> MODE -> TOGGLE BUTTONS
        # These are toggle buttons for selecting the mode of buriedness calculation.
        self._widget_refine_buriedness_mode_buttons = self._create_toggle_buttons_mode(
            observer=self._on_value_change_buriedness_mode
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> BURIEDNESS -> FILTER -> TOGGLE BUTTONS
        # These are toggle buttons for deciding whether the given value range should be included
        # or excluded.
        self._widget_refine_buriedness_filter_buttons = self._create_toggle_buttons_filter(
            observer=self._on_value_change_buriedness_filter
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> BURIEDNESS -> RANGE -> SLIDER
        # This is a range slider for selecting the range of values for vacancy calculation:
        self._widget_refine_buriedness_range_slider = self._create_range_slider(
            observer=self._on_value_change_buriedness_range
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY -> RANGE -> UNIT
        # This is a label that shows the unit of the selected vacancy source values:
        self._widget_refine_buriedness_range_unit = widgets.Label(
            value="",
            layout=widgets.Layout(margin="0px 0px 0px 0px", align="left")
        )

        # CREATE WIDGET FOR: BINDING SITE REFINEMENT -> VACANCY
        # This is a container Box for all vacancy widgets above, plus added labels, making the
        # vacancy control panel:
        buriedness_control_panel = self._create_control_panel(
            controllers=[
                widgets.HBox(self._widget_refine_buriedness_source_buttons),
                widgets.Box(
                    children=[self._widget_refine_buriedness_mode_buttons],
                    layout=widgets.Layout(
                        width='330px',
                        height='',
                        flex_flow='row',
                        display='flex'
                    )
                ),
                self._widget_refine_buriedness_filter_buttons,
                widgets.HBox(
                    [
                        self._widget_refine_buriedness_range_slider,
                        self._widget_refine_buriedness_range_unit
                    ]
                ),
            ],
            header_name="Buriedness",
            header_color="silver"
        )

        # Assemble the refinement panel
        refinement_panel = widgets.HBox([vacancy_control_panel, buriedness_control_panel])

        # Assemble the main panel
        main_panel = widgets.Accordion(
            children=[refinement_panel, widgets.Text(), widgets.Box()],
        )  # layout=Layout(
        # height="200px"))
        main_panel.set_title(0, "Binding Site Refinement")
        main_panel.set_title(1, "Feature Selection")
        main_panel.set_title(2, "Clustering")
        main_panel.selected_index = 0


        self._curr_selection_vacancy_source: np.ndarray = np.zeros(
            shape=self._t2fpharm.fields.fields_count + 1,
            dtype=np.bool_,
        )
        self._curr_selection_vacancy_mode: str = "min"
        self._curr_selection_vacancy_range: Tuple[float, float] = (0, 0)
        self._curr_selection_vacancy_filter: bool = True

        self._curr_values_vacancy: np.ndarray = None
        self._curr_mask_vacancy_lower_bound: np.ndarray = self._initialize_mask_array()
        self._curr_mask_vacancy_upper_bound: np.ndarray = self._initialize_mask_array()
        self._curr_mask_vacancy_range: np.ndarray = self._initialize_mask_array()
        self._curr_mask_vacancy: np.ndarray = self._initialize_mask_array()
        self._curr_mask_buriedness: np.ndarray = self._initialize_mask_array()
        self._curr_mask_refinement: np.ndarray = self._initialize_mask_array()
        self._curr_mask_feature_selection: np.ndarray = self._initialize_mask_array()
        self._curr_mask: np.ndarray = self._initialize_mask_array()

        self._mask_vacancy: np.ndarray
        self._mask_buriedness: np.ndarray
        self._mask_vacancy_and_buriedness: np.ndarray

        self._nglwidget = self._t2fpharm.receptor.create_new_ngl_widget()
        self._update_grid()

        self._debug = widgets.Output()

        display(self._debug, main_panel, self._nglwidget.display(gui=True))
        return

    def _update_grid(self):
        self._curr_mask = self._curr_mask_vacancy
        nglview_api.remove_component_by_name(view=self._nglwidget, name="grid")
        nglview_api.add_spheres(
            view=self._nglwidget,
            coords=self._t2fpharm.fields.grid.coordinates[self._curr_mask[0]],
            name="grid"
        )
        return

    def _initialize_mask_array(self):
        return np.ones(shape=self._t2fpharm.fields.tensor.shape[:-1], dtype=np.bool_)

    def _on_value_change_vacancy_source(self, change: dict):
        source_name: str = change["owner"].description
        selected: bool = change["new"]
        self._curr_selection_vacancy_source[
            -1 if source_name == "Distance" else
            self._t2fpharm.fields.index_field_names(names=[source_name])
        ] = selected
        num_selected_sources = np.count_nonzero(self._curr_selection_vacancy_source)
        self._widget_refine_vacancy_mode_buttons.disabled = num_selected_sources in (0, 1)
        self._widget_refine_vacancy_range_slider.disabled = num_selected_sources == 0
        self._widget_refine_vacancy_filter_buttons.disabled = num_selected_sources == 0
        if num_selected_sources == 0:
            self._curr_mask_vacancy[...] = True
            self._curr_mask_buriedness[...] = True
            self._curr_mask_refinement[...] = True
            self._curr_mask[...] = self._curr_mask_feature_selection

        else:
            self._calculate_vacancy()
        self._update_grid()
        return

    def _on_value_change_vacancy_mode(self, change: dict):
        self._calculate_vacancy()
        self._update_grid()
        return

    #@debounce(1.5)
    def _on_value_change_vacancy_range(self, change: dict):
        # self._calculate_vacancy_stage_range(val_range=change["new"])
        self._calculate_vacancy()
        self._update_grid()
        return

    def _on_value_change_vacancy_filter(self, change: dict):
        change["new"]: str  # Either "Include" or "Exclude"
        self._curr_mask_vacancy = np.logical_not(self._curr_mask_vacancy)
        self._update_grid()
        return

    def _calculate_buriedness_from_begining(self):
        pass

    def _on_value_change_buriedness_source(self, change: dict):
        return

    def _on_value_change_buriedness_mode(self, change: dict):
        return

    def _on_value_change_buriedness_filter(self, change: dict):
        return

    def _on_value_change_buriedness_range(self, change: dict):
        return

    def _calculate_vacancy(self):
        self._curr_mask_vacancy, vacancy_range = vectorized.filter_array(
            array=self._vacancy_data[..., self._curr_selection_vacancy_source],
            minmax_vals=self._widget_refine_vacancy_range_slider.value,
            reduction_op=self._widget_refine_vacancy_mode_buttons.value,
            reduction_axis=-1,
            invert=self._widget_refine_vacancy_filter_buttons.value == "Exclude"
        )
        with self._widget_refine_vacancy_range_slider.hold_trait_notifications():
            self._widget_refine_vacancy_range_slider.min = vacancy_range[0]
            self._widget_refine_vacancy_range_slider.max = vacancy_range[1]
            slider_value = self._widget_refine_vacancy_range_slider.value
            new_slider_value = []
            if vacancy_range[0] <= slider_value[0] <= vacancy_range[1]:
                new_slider_value.append(slider_value[0])
            else:
                new_slider_value.append(vacancy_range[0])
            if vacancy_range[0] <= slider_value[1] <= vacancy_range[1]:
                new_slider_value.append(slider_value[1])
            else:
                new_slider_value.append(vacancy_range[1])
            self._widget_refine_vacancy_range_slider.value = new_slider_value
        return

    @staticmethod
    def _create_toggle_buttons_source(
            labels: Sequence[str],
            tooltips: Union[str, Sequence[str]] = "",
            initial_values: Union[bool, Sequence[bool]] = False,
            button_style: Literal['success', 'info', 'warning', 'danger', ''] = "danger",
            observer: Optional[Callable] = None,
    ):
        if isinstance(tooltips, str):
            tooltips = [tooltips] * len(labels)
        if isinstance(initial_values, bool):
            initial_values = [initial_values] * len(labels)
        toggle_buttons = []
        for label, initial_value, tooltip in zip(labels, initial_values, tooltips):
            toggle_button = widgets.ToggleButton(
                value=initial_value,
                description=label,
                tooltip=tooltip,
                layout=widgets.Layout(width="auto"),
                button_style=button_style
            )
            toggle_buttons.append(toggle_button)
            if observer is not None:
                # Set all buttons to observe the same observer function:
                toggle_button.observe(observer, names="value")
        return toggle_buttons

    @staticmethod
    def _create_toggle_buttons_mode(
            labels: Union[Sequence[str], Sequence[Tuple[str, Any]]] = (
                    ("Min", "min"),
                    ("Avg", "avg"),
                    ("Max", "max"),
                    ("Sum", "sum"),
                    ("All", "all"),
                    ("Any", "any"),
                    ("One", "one")
            ),
            current_value: Any = "min",
            tooltips: Sequence[str] = (
                    'Minimum value of all selected source fields.',
                    'Average value of all selected source fields.',
                    'Maximum value of all selected source fields.',
                    'Sum of all selected source fields.',
                    'All values of selected source fields.',
                    'Any value of selected source fields.',
                    'Only one value of selected source fields.',
            ),
            button_style: Literal['success', 'info', 'warning', 'danger', ''] = "warning",
            disabled: bool = True,
            observer: Optional[Callable] = None,
    ) -> widgets.ToggleButtons:
        toggle_buttons = widgets.ToggleButtons(
            options=labels,
            value=current_value,
            button_style=button_style,
            disabled=disabled,
            tooltips=tooltips,
            style=widgets.ToggleButtonsStyle(button_width="auto"),
        )
        if observer is not None:
            toggle_buttons.observe(observer, names="value")
        return toggle_buttons

    @staticmethod
    def _create_toggle_buttons_filter(
        observer: Optional[Callable] = None,
    ):
        toggle_buttons = widgets.ToggleButtons(
            options=["Include", "Exclude"],
            value="Include",
            disabled=True,
            button_style='info',
            tooltips=['Include the points filtered by current selection.',
                      'Exclude the points filtered by current selection.'],
            style=widgets.ToggleButtonsStyle(button_width="4em"),
        )
        if observer is not None:
            toggle_buttons.observe(observer, names="value")
        return toggle_buttons

    @staticmethod
    def _create_range_slider(
            observer: Optional[Callable] = None,
    ):
        range_slider = widgets.FloatRangeSlider(
            value=[0, 0],
            min=0,
            max=0,
            step=0.01,
            disabled=True,
            continuous_update=False,
            orientation='horizontal',
            readout=True,
            readout_format='.2e',
            layout=widgets.Layout(display="flex", align_items='stretch', width='100%')
        )
        if observer is not None:
            range_slider.observe(observer, names="value")
        return range_slider

    @staticmethod
    def _create_control_panel(
            controllers: Sequence[widgets.Widget],
            header_name: str,
            header_color: str = "lightblue",
            controller_labels: Sequence[str] = ("Source", "Mode", "Filter", "Range"),
    ):
        # Create the header, i.e. title of the widget
        # This is created as a button, because it offers more styling options than text,
        # but it doesn't have any functionality.
        header = widgets.Button(
            description=header_name,
            layout=widgets.Layout(width='auto'),
            style=widgets.ButtonStyle(button_color=header_color, font_weight='bold')
        )
        # Create the left column holding labels of controllers:
        labels_column = widgets.VBox(
            [
                widgets.Label(
                    value=f"{label} :",
                    layout=widgets.Layout(margin="3px 10px 3px 0px", width="50px")
                )
                for label in controller_labels
            ],
            layout=widgets.Layout(width="100px"),
        )
        # Create the right column holding the controllers
        controllers_column = widgets.VBox(
            controllers,
            layout=widgets.Layout(
                display='flex',
                flex_flow='column',
                align_items='stretch',
                width='100%'
            ),
        )

        # Horizontal container for two columns holding labels and controllers:
        controller_box = widgets.HBox([labels_column, controllers_column])

        # Assemble the control panel
        control_panel = widgets.Box(
            layout=widgets.Layout(
                display='flex',
                flex_flow='column',
                align_items='stretch',
                width='50%'
            ),
            children=[header, controller_box],
        )
        return control_panel








# RECEPTOR_FILEPATHS = [Path("/Users/home/Downloads/3w32.pdbqt")]
# GRID_CENTER = (15.91, 32.33, 11.03)
# GRID_SIZE = (20, 20, 20)
# GRID_SPACING = 0.6
# HBD = autodock.AtomType.HD
# HBA = autodock.AtomType.OA
# AL = autodock.AtomType.C
# AR = autodock.AtomType.A
# PROBES = (HBD, HBA, AL, AR)
# VACANCY_MAX_ENERGY = 0.6
# PSP_NUM_DIRS = 7
# PSP_LEN_MAX = 10.0
# PSP_COUNT_MIN = 4
# HBOND_LEN_MAX = 3.0
# t2f = T2FPharm(
#     receptor_filepaths=RECEPTOR_FILEPATHS,
#     ligand_types=PROBES,
#     grid_center=GRID_CENTER,
#     grid_size=GRID_SIZE,
#     grid_spacing=GRID_SPACING
# )
# t2f()



        #
        #
        # if isinstance(receptors, (PathLike, Sequence)):
        #     if isinstance(receptors, PathLike):
        #         path = Path(receptors)
        #         if path.is_file():
        #             if path.suffix.lower() == ".pdbqt":
        #                 pdbqt_filepaths = [path]
        #             else:
        #                 raise ValueError("Receptor's input file's extension should be PDBQT.")
        #         elif path.is_dir():
        #             pdbqt_filepaths = list(path.glob("*.pdbqt"))
        #         else:
        #             raise ValueError(f"No such file or directory: {path}")
        #     else:
        #         pdbqt_filepaths = [
        #             Path(receptor) for receptor in receptors
        #             if Path(receptor).is_file() and Path(receptor).suffix.lower() == ".pdbqt"
        #         ]
        #     if len(pdbqt_filepaths) == 0:
        #         raise ValueError(f"No PDBQT file found.")
        #     if len(pdbqt_filepaths) == 1:
        #         self._receptor = protein.Protein.from_file(path=pdbqt_filepaths[0])
        #         self._receptor_is_static = True
        #         self._receptor_is_trajectory = False
        #     else:
        #         self._receptor_is_static = False
        #         try:
        #             self._receptor = protein.ProteinTrajectory.from_files(paths=pdbqt_filepaths)
        #             self._receptor_is_trajectory: bool = True
        #         except protein.TopologyError:
        #             self._receptor = protein.ProteinFamily.from_files(paths=pdbqt_filepaths)
        #             self._receptor_is_trajectory: bool = False
        # elif isinstance(
        #         receptors,
        #         (protein.Protein, protein.ProteinTrajectory, protein.ProteinFamily)
        # ):
        #     self._receptor = receptors
        #     pdbqt_filepaths = self._receptor.to_file_pdbqt()
        #     if isinstance(receptors, protein.ProteinFamily):
        #         self._receptor_is_static = False
        #         self._receptor_is_trajectory = False
        #     elif isinstance(receptors, protein.ProteinTrajectory):
        #         self._receptor_is_static = False
        #         self._receptor_is_trajectory = True
        #     else:
        #         self._receptor_is_static = True
        #         self._receptor_is_trajectory = False
        #         pdbqt_filepaths = [pdbqt_filepaths]
        # else:
        #     raise ValueError("Type of `receptor` not recognized.")
        #
        #
        # self._output_path = pdbqt_filepaths[0].parent if output_path is None else output_path
        # self._output_path.mkdir(parents=True, exist_ok=True)
        # self._probe_types = np.array(probe_types, dtype=object)
        #
        #
        #
        #
        # self._df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptors[0])["ATOM"]
        #
        # for t, receptor_filepath in enumerate(receptors):
        #     df_pdbqt = pdbqt.parse_pdbqt(filepath_pdbqt=receptor_filepath)["ATOM"]
        #     self._protein_coordinates[t, ...] = df_pdbqt[["x", "y", "z"]].to_numpy()
        #
        # self._grid_dists_protein = np.moveaxis(
        #     distance_matrix(
        #         self._grid_coords.reshape(-1, 3),
        #         self._protein_coordinates.reshape(-1, 3)
        #     ).reshape(*self._grid_shape, self._temporal_len, len(self._df_pdbqt.index)),
        #     source=3,
        #     destination=0,
        # )
        #
        # # Get all 26 half directions, in the order: orthogonal, 3d-diagonal, 2d-diagonal
        # self._direction_vects = np.zeros(shape=(26, 4), dtype=np.byte)
        # self._direction_vects[:, 1:] = np.concatenate(
        #     [spatial.GRID_DIRS_ORTHO, spatial.GRID_DIRS_DIAG_3D, spatial.GRID_DIRS_DIAG_2D]
        # )
        # # Calculate the length of each direction unit vector
        # self._direction_vects_len = np.linalg.norm(self._direction_vects, axis=-1) * self._grid_spacing

