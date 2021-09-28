"""
opencadd.structure.pocket.viewer

Visualize pocket(s).
"""

import logging

from matplotlib import colors
import nglview

from opencadd.io import DataFrame

_logger = logging.getLogger(__name__)


class PocketViewer:
    """
    Defines the NGLview widget to be used to visualize one or more Pocket objects.

    Attributes
    ----------
    viewer : nglview.widget.NGLWidget
        NGLview widget. Call this if you want to view pockets in e.g. Jupyter Notebooks.
    pockets_residue_ngl_ixs : dict of (int or str: list of int)
        For each structure (key), list of pocket NGLview indices.
    structure_names : list of (int or str)
        Names for all structures/pockets in the viewer.
    _residue_ids_to_ngl_ixs = dict of pandas.DataFrames
        For each structure (keys), mapping of all residue PDB IDs and names to the NGLview residue
        indices (values).
    _component_counter : int
        Next component index to be filled with a component.
    _components_structures : dict of int
        For each structure (keys), the NGLview component index.
        - Structure key names: <pocket name> (taken from Pocket attribute "name")
    _components_pocket_center : dict of int
        For each structure (keys), the NGLview component index.
        - Structure key names: <pocket name> (taken from Pocket attribute "name")
    _components_subpockets : dict of dict of int
        For each structure (keys - level 1) and for each subpocket sphere (keys - level 2),
        the NGLview component index.
        - Structure key names: <pocket name> (taken from Pocket attribute "name")
        - Subpocket key names: <subpocket name> (taken from Subpocket attribute "name")
    _components_anchor_residues :
        For each structure (keys - level 1) and for each anchor residue sphere (keys - level 2),
        the NGLview component index.
        - Structure key names: <pocket name> (taken from Pocket attribute "name")
        - Anchor residue key names: <subpocket name>_<anchor residue PDB ID>
          (taken from Subpocket attribute "name" and AnchorResidue attribute "residue_id")
    """

    def __init__(self):

        self.viewer = nglview.NGLWidget()
        self.viewer._remote_call("setSize", target="Widget", args=["1000px", "600px"])
        self.pockets_residue_ngl_ixs = {}
        self.structure_names = []
        self._residue_ids_to_ngl_ixs = {}
        self._component_counter = 0
        self._components_structures = {}
        self._components_pocket_center = {}
        self._components_subpockets = {}
        self._components_anchor_residues = {}

    def add_pocket(
        self,
        pocket,
        ligand_expo_id=None,
        show_pocket_center=True,
        show_subpockets=True,
        show_anchor_residues=True,
        show_regions=True,
        sphere_opacity=0.7,
        sphere_color_pocket_center="blue",
        show_only_pocket_residues=False,
    ):
        """
        Visualize the pocket (subpockets, regions, and anchor residues).

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object to be visualized.
        ligand_expo_id : str or None
            Expo ID for co-crystallized ligand as defined in PDB.
            If None (default), ligand is not shown.
        show_pocket_center : bool
            Show the pocket center as sphere (default) or not.
        show_subpockets : bool
            Show the subpockets as spheres (default) or not.
        show_anchor_residues : bool
            Show the anchor residues as spheres (default) or not.
        show_regions : bool
            Color residues by regions (default) or not.
        sphere_opacity : float
            Sphere opacity (default 0.7). Applied to all spheres.
        sphere_color_pocket_center : str
            Pocket center sphere color (matplotlib name). Default is blue.

        Returns
        -------
        nglview.widget.NGLWidget
            Pocket visualization.
        """

        # Get residue ID > nglview index mapping based on file format
        # (this is important to know how nglview will index residues)
        self._map_residue_ids_names_nglixs(pocket)
        self.pockets_residue_ngl_ixs[pocket.name] = (
            pocket.residues.dropna()
            .merge(self._residue_ids_to_ngl_ixs[pocket.name], how="left", on="residue.id")[
                "residue.ngl_ix"
            ]
            .to_list()
        )

        # Load structure from text in nglview
        self._add_structure(pocket, ligand_expo_id)

        # Show regions
        if show_regions:
            self._add_regions(pocket, show_only_pocket_residues)

        # Show pocket center
        if show_pocket_center:
            self._add_pocket_center(pocket, sphere_opacity, sphere_color_pocket_center)

        # Show subpockets
        if show_subpockets:
            self._add_subpockets(pocket, sphere_opacity)

        # Show anchor points
        if show_anchor_residues:
            self._add_anchor_residues(pocket, sphere_opacity)

        # Show
        return self.viewer

    def hide(self, structure_name):
        """
        Hide a structure and all related components (such as spheres).
        The components are only hidden, not removed. They can be shown again using ".show()".

        Parameters
        ----------
        structure_name : str
            Structure/pocket name (as defined in the Pocket object ("Pocket.name")).
        """

        components_to_be_hidden = self._components_by_structure_name(structure_name)
        self.viewer.hide(components_to_be_hidden)

    def show_only(self, structure_name):
        """
        Show a structure and all related components (such as spheres).

        Parameters
        ----------
        structure_name : str
            Structure/pocket name (as defined in the Pocket object ("Pocket.name")).
        """

        components_to_be_shown = self._components_by_structure_name(structure_name)
        self.viewer.show_only(components_to_be_shown)

    def show_all(self):
        """
        Show all components.
        """

        self.viewer.show()

    def _map_residue_ids_names_nglixs(self, pocket):
        """
        Map residue IDs and names to nglview indices depending on file format.
        In case of mol2 files, nglview will use indices starting from 1.
        In case of pdb files, nglview will use the residue IDs as indices.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.

        Returns
        -------
        pandas.Series
            Residue IDs (index) and residue nglview indices (values).
        """

        # Get atom data
        # Cast residue IDs to integer - drop atoms where this is not possible!
        dataframe = DataFrame.from_text(pocket._text, pocket._extension)
        drop_ixs = []
        for index, residue_id in dataframe["residue.id"].items():
            try:
                residue_id = int(residue_id)
            except (TypeError, ValueError):
                drop_ixs.append(index)
        dataframe.drop(drop_ixs, inplace=True)
        dataframe = dataframe.astype({"residue.id": "int32"})

        # Get all residue names and IDs (full structure!!)
        residue_id2ix = dataframe[["residue.name", "residue.id"]].drop_duplicates()

        if pocket._extension == "mol2":

            # Map residue names to nglview index (starting from 1)
            residue_id2ix["residue.ngl_ix"] = [str(i) for i in range(1, len(residue_id2ix) + 1)]

        else:

            # In this case, residue ID and nglview index are the same
            residue_id2ix["residue.ngl_ix"] = [str(i) for i in residue_id2ix["residue.id"]]

        self._residue_ids_to_ngl_ixs[pocket.name] = residue_id2ix

    def _add_structure(self, pocket, ligand_expo_id=None):
        """
        Add protein structure as defined in pocket object ("_text" attribute) in cartoon
        representation. Optionally, show ligand in ball+stick representation.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.
        ligand_expo_id : str or None
            Expo ID for co-crystallized ligand as defined in PDB.
            If None (default), ligand is not shown.
        """

        # Set structure name
        if pocket.name in self.structure_names:
            raise ValueError(
                f"Structure/pocket name ({pocket.name}) already taken. "
                f"FYI: The following structures are already in your view: {self.structure_names}"
            )
        self.structure_names.append(pocket.name)
        # Add structure
        component = nglview.TextStructure(pocket._text, ext=pocket._extension)
        self.viewer.add_component(component)
        # Save the structure's NGLview component
        self._components_structures[pocket.name] = self._component_counter
        self._component_counter += 1
        # Clear the structure's representation
        self.viewer.clear_representations(component=self._components_structures[pocket.name])
        # Add protein cartoon representation
        self.viewer.add_representation(
            "cartoon",
            selection="protein",
            component=self._components_structures[pocket.name],
            color="grey",
        )
        # If specified, add ligand ball+stick representation
        if ligand_expo_id:
            self._add_ligand(pocket, ligand_expo_id)

    def _add_ligand(self, pocket, ligand_expo_id):
        """
        Show ligand in ball+stick representation.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.
        ligand_expo_id : str
            Expo ID for co-crystallized ligand as defined in PDB.
        """

        # Cast residue mapping to Series (residue name as index, NGL index as values)
        residue_name2ix = self._residue_ids_to_ngl_ixs[pocket.name]
        residue_name2ix = residue_name2ix.set_index("residue.name")["residue.ngl_ix"]

        try:
            residue_ngl_ix = residue_name2ix.loc[ligand_expo_id]
            self.viewer.add_representation(
                "ball+stick",
                selection=residue_ngl_ix,
                component=self._components_structures[pocket.name],
            )
        except KeyError:
            if pocket._extension != "mol2":
                self.viewer.add_representation(
                    "ball+stick",
                    selection="hetero and not water and not ions",
                    component=self._components_structures[pocket.name],
                )
            _logger.warning(
                f"Ligand {ligand_expo_id} could not be found. "
                f"If you are using the mol2 format, the ligand name might not be read correctly. "
                f"We will try to infer and show the ligand based on the following selection: "
                f"hetero and not water and not ions"
            )

    def _add_regions(self, pocket, show_only_pocket_residues=False):
        """
        Color residues by regions.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.
        """

        # Do nothing if pocket has no regions
        if pocket.regions is None:
            _logger.info(f"Pocket {pocket.name} has no regions.")
            return None

        # Cast residue mapping to Series (residue IDs as index, NGL index as values)
        residue_id2ix = self._residue_ids_to_ngl_ixs[pocket.name]
        residue_id2ix = residue_id2ix.set_index("residue.id")["residue.ngl_ix"]

        regions = pocket.regions.dropna(axis=0, subset=["residue.id"])
        scheme_regions_list = []
        for _, region in regions.iterrows():
            color = region["region.color"]
            residue_id = region["residue.id"]
            residue_ngl_ix = residue_id2ix.loc[residue_id]
            scheme_regions_list.append([color, residue_ngl_ix])
        scheme_regions = nglview.color._ColorScheme(scheme_regions_list, label="scheme_regions")
        if show_only_pocket_residues:
            selection = " or ".join(self.pockets_residue_ngl_ixs[pocket.name])
        else:
            selection = "protein"
        self.viewer.clear_representations(self._components_structures[pocket.name])
        self.viewer.add_representation(
            "cartoon",
            selection=selection,
            component=self._components_structures[pocket.name],
            color=scheme_regions,
        )

    def _add_pocket_center(self, pocket, sphere_opacity=0.7, sphere_color="blue"):
        """
        Add the pocket center as sphere.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.
        sphere_opacity : float
            Sphere opacity (default 0.7).
        sphere_color : str
            Sphere color (matplotlib name). Default ist blue.
        """

        # Do nothing if pocket has no center
        if pocket.center is None:
            _logger.info(f"Pocket {pocket.name} has no pocket center.")
            return None

        color_rgb = colors.to_rgb(sphere_color)
        # Cast numpy.float to float (otherwise TypeError!)
        center = [float(i) for i in pocket.center]
        self.viewer.shape.add_sphere(center, color_rgb, 2, f"center: {pocket.name}")
        # Save NGLview component
        self._components_pocket_center[pocket.name] = self._component_counter
        self._component_counter += 1
        # Set sphere opacity
        self.viewer.update_representation(
            component=self._components_pocket_center[pocket.name],
            repre_index=0,
            opacity=sphere_opacity,
        )

    def _add_subpockets(self, pocket, sphere_opacity=0.7):
        """
        Add the pocket's subpockets as spheres.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.
        sphere_opacity : float
            Sphere opacity (default 0.7).
        """

        # Do nothing if pocket has no subpockets
        if pocket.subpockets is None:
            _logger.info(f"Pocket {pocket.name} has no subpockets.")
            return None

        self._components_subpockets[pocket.name] = {}
        for subpocket in pocket._subpockets:
            if subpocket.center is not None:
                # Cast numpy.float to float (otherwise TypeError!)
                center = [float(i) for i in subpocket.center]
                name = subpocket.name
                color_rgb = colors.to_rgb(subpocket.color)
                self.viewer.shape.add_sphere(center, color_rgb, 2, f"{name}: {pocket.name}")
                # Save NGLview component
                self._components_subpockets[pocket.name][subpocket.name] = self._component_counter
                self._component_counter += 1
                # Set sphere opacity
                self.viewer.update_representation(
                    component=self._components_subpockets[pocket.name][subpocket.name],
                    repre_index=0,
                    opacity=sphere_opacity,
                )

    def _add_anchor_residues(self, pocket, sphere_opacity=0.7):
        """
        Add the subpocket's anchor residues as spheres.

        Parameters
        ----------
        pocket : opencadd.structure.pocket.Pocket
            Pocket object.
        sphere_opacity : float
            Sphere opacity (default 0.7).
        """

        # Do nothing if pocket has no anchor residues
        if pocket.anchor_residues is None:
            _logger.info(f"Pocket {pocket.name} has no subpockets.")
            return None

        self._components_anchor_residues[pocket.name] = {}
        for _, anchor_residue in pocket.anchor_residues.iterrows():
            if anchor_residue["anchor_residue.center"] is not None:
                # Cast numpy.float to float (otherwise TypeError!)
                center = [float(i) for i in anchor_residue["anchor_residue.center"]]
                color_rgb = colors.to_rgb(anchor_residue["anchor_residue.color"])
                self.viewer.shape.add_sphere(center, color_rgb, 0.5)
                # Save NGLview component
                self._components_anchor_residues[pocket.name][
                    f"{anchor_residue['subpocket.name']}_{anchor_residue['anchor_residue.id']}"
                ] = self._component_counter
                self._component_counter += 1
                # Set sphere opacity
                self.viewer.update_representation(
                    component=self._components_anchor_residues[pocket.name][
                        f"{anchor_residue['subpocket.name']}_{anchor_residue['anchor_residue.id']}"
                    ],
                    repre_index=0,
                    opacity=sphere_opacity,
                )

    def _components_by_structure_name(self, structure_name):
        """
        Get all components that are connected to a given structure.

        Parameters
        ----------
        structure_name : str
            Structure/pocket name.
        """

        if structure_name not in self.structure_names:
            raise ValueError(
                f"Structure/pocket name ({structure_name}) does not exist. "
                f"FYI: The following structures are in your view: {self.structure_names}"
            )

        components = []
        components.append(self._components_structures[structure_name])
        components.append(self._components_pocket_center[structure_name])
        if len(self._components_subpockets) > 0:
            components.extend(list(self._components_subpockets[structure_name].values()))
        if len(self._components_anchor_residues) > 0:
            components.extend(list(self._components_anchor_residues[structure_name].values()))

        return components
