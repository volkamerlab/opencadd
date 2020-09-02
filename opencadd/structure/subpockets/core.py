"""
core.py

Defines core classes and functions.
"""

from matplotlib import colors
import nglview
import numpy as np
import pandas as pd


class Base:
    """
    Basic methods needed for child classes in this module.
    """

    def _format_color(self, color):
        """
        Handle input color name: Return color name and RGB value.

        Parameters
        ----------
        color : str
            Color name.

        Returns
        -------
        tuple of (str, list of float)
            Color name, color RGB value.
        """

        if isinstance(color, str):
            color_name = color
            color_rgb = colors.to_rgb(color_name)

        return color_name, color_rgb

    def _format_residue_pdb_ids_and_labels(self, residue_pdb_ids, residue_labels):
        """
        Handle input residue PDB IDs and labels: Must be of same length, cast values to string.

        Parameters
        ----------

        Returns
        -------
        tuple of list of str
            Residue PDB IDs, residue labels.
        """

        # If no anchor residue labels are given, create list of None
        # (list length = number of anchor residue PDB IDs)
        if not residue_labels:
            residue_labels = [None] * len(residue_pdb_ids)

        # Check if PDB IDs and labels are of same length, if not raise error
        if len(residue_pdb_ids) != len(residue_labels):
            raise ValueError(f"Number of residue PDB IDs and labels must be of same length.")

        # Cast residue PDB IDs and labels to strings (except None)
        residue_pdb_ids = [
            str(residue_pdb_id) if residue_pdb_id else residue_pdb_id
            for residue_pdb_id in residue_pdb_ids
        ]
        residue_labels = [
            str(residue_label) if residue_label else residue_label
            for residue_label in residue_labels
        ]

        return residue_pdb_ids, residue_labels


class Pocket(Base):
    """
    Class defining a pocket with 
    - structural protein data, 
    - subpockets (to be shown as spheres) and 
    - regions (to be highlighted).

    Attributes
    ----------
    data : pandas.DataFrame
        Structural protein data with the following mandatory columns:
        "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z".
    name : str
        Name of protein.
    residue_pdb_ids : list of str
        Pocket residue PDB IDs.
    residue_labels : list of str
        Pocket residue labels.
    centroid : numpy.array
        Pocket centroid based on all pocket residues' CA atoms.
    _subpockets : list of Subpocket
        List of user-defined subpockets.
    _region : list of Region
        List of user-defined regions.
    """

    def __init__(self, data, name, residue_pdb_ids, residue_labels):

        self.data = data
        self.name = name
        residue_pdb_ids, residue_labels = self._format_residue_pdb_ids_and_labels(
            residue_pdb_ids, residue_labels
        )
        self.residue_pdb_ids = residue_pdb_ids
        self.residue_labels = residue_labels
        self.centroid = self._centroid()
        self._subpockets = []
        self._regions = []

    @property
    def residues(self):
        """
        Return pocket residues (PDB ID and labels).
        """

        residues = {"residue.pdb_id": self.residue_pdb_ids, "residue.label": self.residue_labels}
        residues = pd.DataFrame(residues)
        return residues.reset_index(drop=True)

    @property
    def subpockets(self):
        """
        Return subpockets data as DataFrame: 
        Name, color and subpocket center (columns) for all subpockets (rows).
        """

        if self._subpockets == []:
            return None

        subpockets = pd.DataFrame(
            {
                "subpocket.name": [subpocket.name for subpocket in self._subpockets],
                "subpocket.color_name": [subpocket.color_name for subpocket in self._subpockets],
                "subpocket.color_rgb": [subpocket.color_rgb for subpocket in self._subpockets],
                "subpocket.center": [subpocket.center for subpocket in self._subpockets],
            }
        )

        return subpockets.reset_index(drop=True)

    @property
    def regions(self):
        """
        Return region data as DataFrame:
        Name, color, involved residue PDB IDs and labels (columns) for all regions.
        """
        if self._regions == []:
            return None

        regions = []

        for region in self._regions:

            n_residues = len(region.residue_pdb_ids)

            region = pd.DataFrame(
                {
                    "region.name": [region.name] * n_residues,
                    "region.color_name": [region.color_name] * n_residues,
                    "region.color_rgb": [region.color_rgb] * n_residues,
                    "residue.pdb_id": region.residue_pdb_ids,
                    "residue.label": region.residue_labels,
                }
            )

            regions.append(region)

        regions = pd.concat(regions)

        return regions.reset_index(drop=True)

    @property
    def anchor_residues(self):
        """
        Return anchor residue data as DataFrame:
        - Subpocket name and color
        - Anchor residue PDB IDs (user-defined input IDs or alternative 
          IDs if input was not available)
        - Anchor residue labels
        - Ahe anchor residue centers (coordinates)
        """

        if self._subpockets == []:
            return None

        anchor_residues = [subpocket.anchor_residues for subpocket in self._subpockets]
        anchor_residues = pd.concat(anchor_residues)

        return anchor_residues.reset_index(drop=True)

    def clear_subpockets(self):
        """
        Clear subpockets, i.e. remove all defined subpockets.
        """

        self._subpockets = []

    def clear_regions(self):
        """
        Clear regions, i.e. remove all defined regions.
        """

        self._regions = []

    def _centroid(self):
        """
        Add centroid of all input residues' CA atoms.

        Returns
        ----------
        numpy.array
            Pocket centroid (coordinates).
        """

        dataframe = self.data

        atoms = dataframe[
            (dataframe["residue.pdb_id"].isin(self.residue_pdb_ids))
            & (dataframe["atom.name"] == "CA")
        ]

        print(f"The pocket centroid is calculated based on {len(atoms)} CA atoms.")

        centroid = atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()

        return centroid

    def add_subpocket(
        self, name, color, anchor_residue_pdb_ids, anchor_residue_labels=None,
    ):
        """
        Add subpocket based on given anchor residue PDB IDs.

        Parameters
        ----------
        name : str
            Subpocket name.
        color : str
            Subpocket color.
        anchor_residue_pdb_ids : list of (int, str)
            List of anchor residue PDB IDs.
        anchor_residue_labels : list of (int, str) or None
            List of anchor residue labels. Must be of same length as anchor_residue_pdb_ids.
        """

        subpocket = Subpocket()
        subpocket.from_dataframe(
            self.data, name, color, anchor_residue_pdb_ids, anchor_residue_labels
        )
        self._subpockets.append(subpocket)

    def add_region(self, name, color, residue_pdb_ids, residue_labels=None):
        """
        Add region based on given input residue PDB IDs.

        Parameters
        ----------
        name : str
            Region name.
        color : str
            Region color.
        residue_pdb_ids : list of (int, str)
            List of residue PDB IDs defining the region.
        residue_labels : list of (int, str) or None
            List of residue labels. Must be of same length as residue_pdb_ids.
        """

        region = Region()
        region.from_dataframe(self.data, name, color, residue_pdb_ids, residue_labels)
        self._regions.append(region)

    def visualize(self, filepath):
        """
        Visualize the pocket (subpockets, regions, and anchor residues).

        Parameters
        ----------
        filepath : pathlib.Path or str
            Path to structure file.
        
        Returns
        -------
        nglview.widget.NGLWidget
            Pocket visualization.
        """

        filepath = str(filepath)

        # Load structure from file in nglview
        view = nglview.show_file(filepath)
        view._remote_call("setSize", target="Widget", args=["1000px", "600px"])
        view.clear()

        # Get file format (this is important to know how nglview will index residues)
        file_format = filepath.split(".")[-1]

        # Get PDB ID to nglview index mapping
        residue_id2ix = self._map_residue_id2ix(file_format)

        # Show regions
        scheme_regions_list = []
        for index, region in self.regions.iterrows():
            color_name = region["region.color_name"]
            residue_pdb_id = region["residue.pdb_id"]
            residue_ngl_ix = residue_id2ix.loc[residue_pdb_id]
            scheme_regions_list.append([color_name, residue_ngl_ix])

        scheme_regions = nglview.color._ColorScheme(scheme_regions_list, label="scheme_regions")
        view.add_representation("cartoon", selection="protein", color=scheme_regions)

        # Show pocket centroids
        view.shape.add_sphere(list(self.centroid), [0, 0, 1], 2, "centroid")

        # Show subpockets
        for index, subpocket in self.subpockets.iterrows():
            center = list(subpocket["subpocket.center"])
            color_rgb = subpocket["subpocket.color_rgb"]
            name = subpocket["subpocket.name"]
            view.shape.add_sphere(center, color_rgb, 2, name)

        # Show anchor points
        for index, anchor_residue in self.anchor_residues.iterrows():
            center = list(anchor_residue["anchor_residue.center"])
            color_rgb = anchor_residue["subpocket.color_rgb"]
            view.shape.add_sphere(center, color_rgb, 0.5)

        # Show
        return view

    def _map_residue_id2ix(self, file_format):
        """
        Map residue PDB IDs to nglview indices depending on file format.
        In case of mol2 files, nglview will use indices starting from 1.
        In case of pdb files, nglview will use the residue PDB IDs as indices.

        Parameters
        ----------
        filepath : pathlib.Path or str
            Path to structure file.

        Returns
        -------
        pandas.Series
            Residue PDB IDs (index) and residue nglview indices (values).
        """

        # Get all residue names
        residue_id2ix = self.data[["residue.subst_name", "residue.pdb_id"]].drop_duplicates()

        if file_format == "mol2":

            # Map residue names to nglview index (starting from 1)
            residue_id2ix["residue.ngl_ix"] = [str(i) for i in range(1, len(residue_id2ix) + 1)]
            # Cast to Series (PDB IDs as index, NGL index as values)
            residue_id2ix.set_index("residue.pdb_id", inplace=True)
            residue_id2ix = residue_id2ix["residue.ngl_ix"]

        else:

            # In this case, PDB ID and nglview index are the same
            residue_id2ix["residue.ngl_ix"] = residue_id2ix["residue.pdb_id"]
            residue_id2ix.index = residue_id2ix["residue.pdb_id"]
            residue_id2ix = residue_id2ix["residue.ngl_ix"]

        return residue_id2ix


class Subpocket(Base):
    """
    Class defining a subpocket.

    Attributes
    ----------
    name : str
        Subpocket name.
    color_name : str
        Region color name.
    color_rgb : list of float
        Region color RGB value.
    center : np.array
        Coordinates (x, y, z) of the subpocket center, 
        i.e. the centroid of all anchor residues' CA atoms.
    _anchor_residues : list of Residue
        List of anchor residues.
    """

    def __init__(self):

        self.name = None
        self.color_name = None
        self.color_rgb = None
        self.center = None
        self._anchor_residues = None

    @property
    def anchor_residues(self):
        """
        Return anchor residue data as DataFrame:
        - Subpocket name and color
        - Anchor residue PDB IDs (user-defined input IDs or alternative 
          IDs if input was not available)
        - Anchor residue labels
        - Ahe anchor residue centers (coordinates)
        """

        anchor_residues_dict = {
            "subpocket.color_name": [residue.color_name for residue in self._anchor_residues],
            "subpocket.color_rgb": [residue.color_rgb for residue in self._anchor_residues],
            "anchor_residue.pdb_id": [residue.pdb_id for residue in self._anchor_residues],
            "anchor_residue.pdb_id_alternative": [
                residue.pdb_id_alternative for residue in self._anchor_residues
            ],
            "anchor_residue.label": [residue.label for residue in self._anchor_residues],
            "anchor_residue.center": [residue.center for residue in self._anchor_residues],
        }
        anchor_residues_df = pd.DataFrame(anchor_residues_dict)
        anchor_residues_df.insert(0, "subpocket.name", self.name)

        return anchor_residues_df

    def from_dataframe(
        self, dataframe, name, color, anchor_residue_pdb_ids, anchor_residue_labels
    ):
        """
        Set subpocket properties.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural protein data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z".
        name : str
            Subpocket name.
        color : str
            Subpocket color (by name or RGB values).
        anchor_residue_pdb_ids : list of (int, str)
            List of anchor residue PDB IDs.
        anchor_residue_labels : list of (int, str) or None
            List of anchor residue labels. Must be of same length as residue_pdb_ids.
        """

        # Format residue PDB IDs and labels
        anchor_residue_pdb_ids, anchor_residue_labels = self._format_residue_pdb_ids_and_labels(
            anchor_residue_pdb_ids, anchor_residue_labels
        )

        anchor_residues = []

        for residue_pdb_id, residue_label in zip(anchor_residue_pdb_ids, anchor_residue_labels):

            residue = AnchorResidue()
            residue.from_dataframe(dataframe, residue_pdb_id, residue_label, color)
            anchor_residues.append(residue)

        self._from_anchor_residues(name, color, anchor_residues)

    def _from_anchor_residues(self, name, color, anchor_residues):
        """
        Set subpocket from given anchor residues.

        Parameters
        ----------
        name : str
            Subpocket name.
        color : str
            Subpocket color.
        anchor_residues : list of Residue
            List of anchor residues.
        """

        self.name = name
        self._anchor_residues = anchor_residues

        # Format and set color name and RGB value
        color_name, color_rgb = self._format_color(color)
        self.color_name = color_name
        self.color_rgb = color_rgb

        # Calculate subpocket center
        self.center = self._centroid()

    def _centroid(self):
        """
        Calculate the centroid of given input anchor residue centers.

        Returns
        -------
        np.array
            Subpocket center, i.e. the centroid of all anchor residue centers. 
            None if anchor residues are missing.
        """

        anchor_residue_centers = [
            anchor_residue.center for anchor_residue in self._anchor_residues
        ]
        # Are there empty anchor residue centers?
        anchor_residue_centers_none = [
            center for center in anchor_residue_centers if center is None
        ]
        # If so, do not return a subpocket center.
        if len(anchor_residue_centers_none) != 0:
            return None
        # Else, calculate the centroid of all given anchor residue centers.
        else:
            subpocket_center = np.mean(anchor_residue_centers, axis=0)
            return subpocket_center


class Region(Base):
    """
    Class defining a region.

    Attributes
    ----------
    name : str
        Region name.
    color_name : str
        Region color name.
    color_rgb : list of float
        Region color RGB value.
    residue_pdb_ids : list of (int, str)
        List of residue PDB IDs defining the region.
    residue_labels : list of (int, str)
        List of residue labels.
    """

    def __init__(self):

        self.name = None
        self.color_name = None
        self.color_rgb = None
        self.residue_pdb_ids = None
        self.residue_labels = None

    def from_dataframe(self, dataframe, name, color, residue_pdb_ids, residue_labels=None):
        """
        Set region properties.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z"
        name : str
            Region name.
        color : str
            Region color.
        residue_pdb_ids : list of (int, str)
            List of residue PDB IDs defining the region.
        residue_labels : list of (int, str) or None
            List of residue labels. Must be of same length as residue_pdb_ids.
        """

        self.name = name

        # Format and set color name and RGB value
        color_name, color_rgb = self._format_color(color)
        self.color_name = color_name
        self.color_rgb = color_rgb

        # Format residue PDB IDs and labels
        residue_pdb_ids, residue_labels = self._format_residue_pdb_ids_and_labels(
            residue_pdb_ids, residue_labels
        )

        # Add residue labels to dataframe
        residue_labels_df = pd.DataFrame(
            {"residue.pdb_id": residue_pdb_ids, "residue.label": residue_labels}
        )
        dataframe = dataframe.merge(residue_labels_df, on="residue.pdb_id", how="left")

        # Keep only existing residue PDB IDs
        residues = dataframe[["residue.pdb_id", "residue.label"]].drop_duplicates()
        residues.reset_index(drop=True, inplace=True)
        residues = residues[residues["residue.pdb_id"].isin(residue_pdb_ids)]
        self.residue_pdb_ids = residues["residue.pdb_id"].to_list()
        self.residue_labels = residues["residue.label"].to_list()


class AnchorResidue(Base):
    """
    Class defining an anchor residue.

    Attributes
    ----------
    pdb_id : str
        Residue PDB ID.
    pdb_id_alternative : list of str
        Alternative residue PDB ID(s) in case input PDB ID not available.
    label : str
        Residue label, e.g. some non-PDB ID.
    color_name : str
        Residue color name.
    color_rgb : list of float
        Residue color RGB value.
    center : numpy.array
        Coordinates (x, y, z) of the residue center.
    """

    def __init__(self):

        self.pdb_id = None
        self.pdb_id_alternative = None
        self.label = None
        self.color_name = None
        self.color_rgb = None
        self.center = None

    def from_dataframe(self, dataframe, residue_pdb_id, residue_label=None, color="green"):
        """
        Set residue properties.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z"
        residue_pdb_id : str
            Residue PDB ID.
        residue_label : str
            Residue label, e.g. some non-PDB ID (default None).
        color : str
            Residue color (default green).
        """

        # Set class attributes
        if isinstance(residue_pdb_id, int):
            residue_pdb_id = str(residue_pdb_id)

        self.pdb_id = residue_pdb_id
        self.label = residue_label

        # Set color name and RGB value
        color_name, color_rgb = self._format_color(color)
        self.color_name = color_name
        self.color_rgb = color_rgb

        # Select atom from residue PDB ID and atom name
        atom = self._ca_atom(dataframe)

        # If residue PDB ID exists, get atom coordinates.
        if atom is not None:
            self.center = atom[["atom.x", "atom.y", "atom.z"]].squeeze().to_numpy()

        # If not, get atom coordinates for residue before/after.
        else:
            atoms = self._ca_atom_before_after(dataframe)
            if atoms is not None:
                self.pdb_id_alternative = atoms["residue.pdb_id"].to_list()
                self.center = atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()
            else:
                self.center = None

    def _ca_atom(self, dataframe):
        """
        Select a CA atom based on a residue PBD ID.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z"

        Returns
        -------
        pandas.DataFrame
            Atom data if DataFrame length is 1, None if length is 0. 
        
        Raises
        ------
        ValueError
            If returned number of atoms is larger than 1.
        """

        atom = dataframe[
            (dataframe["residue.pdb_id"] == self.pdb_id) & (dataframe["atom.name"] == "CA")
        ]

        if len(atom) == 1:
            return atom
        elif len(atom) == 0:
            return None
        else:
            raise ValueError(
                f"Unambiguous atom selection. {len(atom)} atoms found instead of 0 or 1."
            )

    def _ca_atom_before_after(self, dataframe):
        """
        Select CA atoms from residues before and after a given residue PBD ID.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z"

        Returns
        -------
        pandas.DataFrame
            Atoms data if DataFrame length is 1 or 2, None if length is 0. 
        
        Raises
        ------
        ValueError
            If returned number of atoms is larger than 2.
        """

        residue_pdb_id_before = str(int(self.pdb_id) - 1)
        residue_pdb_id_after = str(int(self.pdb_id) + 1)

        atoms = dataframe[
            (dataframe["residue.pdb_id"].isin([residue_pdb_id_before, residue_pdb_id_after]))
            & (dataframe["atom.name"] == "CA")
        ]

        if len(atoms) == 2:
            return atoms
        elif len(atoms) == 1:
            return atoms
        elif len(atoms) == 0:
            return None
        else:
            raise ValueError(
                f"Unambiguous atom selection. {len(atom)} atoms found instead of 0 or 1."
            )
