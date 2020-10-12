"""
opencadd.structure.pocket.pocket

Defines pockets.
"""

from matplotlib import colors
import nglview
import pandas as pd

from .core import Base
from .region import Region
from .subpocket import Subpocket
from opencadd.io import DataFrame


class Pocket(Base):
    """
    Class defining a pocket with
    - structural protein data,
    - subpockets (to be shown as spheres) and
    - regions (to be highlighted).

    Attributes
    ----------
    filepath : str or pathlib.Path
        File path to structural protein data.
    data : pandas.DataFrame
        Structural protein data with the following mandatory columns:
        "residue.id", "atom.name", "atom.x", "atom.y", "atom.z".
    name : str
        Name of protein.
    residue_ids : list of str
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

    def __init__(self):

        self.filepath = None
        self.data = None
        self.name = None
        self.residue_ids = None
        self.residue_labels = None
        self.centroid = None
        self._subpockets = []
        self._regions = []

    @classmethod
    def from_file(cls, filepath, residue_ids, name="", residue_labels=None):
        """
        Initialize Pocket object from structure protein file.

        Attributes
        ----------
        filepath : str or pathlib.Path
            File path to structural protein data.
        residue_ids : list of str
            Pocket residue PDB IDs.
        name : str
            Name of protein (default: empty string).
        residue_labels : None or list of str
            Pocket residue labels. Set to None by default.

        Returns
        -------
        Pocket
            Pocket object.
        """

        pocket = cls()

        pocket.filepath = filepath
        pocket.data = DataFrame.from_file(filepath)
        pocket.name = name
        residue_ids, residue_labels = pocket._format_residue_ids_and_labels(
            residue_ids, residue_labels
        )
        pocket.residue_ids = residue_ids
        pocket.residue_labels = residue_labels
        pocket.centroid = pocket._centroid()

        return pocket

    @property
    def residues(self):
        """
        All pocket's residues.

        Returns
        -------
        pandas.DataFrame
            Residue ID and residue label (columns) for all pocket residues (rows).
        """

        residues = {"residue.id": self.residue_ids, "residue.label": self.residue_labels}
        residues = pd.DataFrame(residues)
        return residues.reset_index(drop=True)

    @property
    def subpockets(self):
        """
        All pocket's subpockets.

        Returns
        -------
        pandas.DataFrame
            Name, color and subpocket center (columns) for all subpockets (rows).
        """

        if not self._subpockets:
            return None

        subpockets = pd.DataFrame(
            {
                "subpocket.name": [subpocket.name for subpocket in self._subpockets],
                "subpocket.color": [subpocket.color for subpocket in self._subpockets],
                "subpocket.center": [subpocket.center for subpocket in self._subpockets],
            }
        )

        return subpockets.reset_index(drop=True)

    @property
    def regions(self):
        """
        All pocket's regions.

        Returns
        -------
        pandas.DataFrame
            Name, color, involved residue PDB IDs and labels (columns) for all regions.
        """
        if self._regions == []:
            return None

        regions = []

        for region in self._regions:

            n_residues = len(region.residue_ids)

            region = pd.DataFrame(
                {
                    "region.name": [region.name] * n_residues,
                    "region.color": [region.color] * n_residues,
                    "residue.id": region.residue_ids,
                    "residue.label": region.residue_labels,
                }
            )

            regions.append(region)

        regions = pd.concat(regions)

        return regions.reset_index(drop=True)

    @property
    def anchor_residues(self):
        """
        All pocket's anchor residues.

        Returns
        -------
        pandas.DataFrame
            Anchor residues (rows) with the following columns:
            - Subpocket name and color
            - Anchor residue PDB IDs (user-defined input IDs or alternative
            IDs if input was not available)
            - Anchor residue labels
            - The anchor residue centers (coordinates)
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
            (dataframe["residue.id"].isin(self.residue_ids)) & (dataframe["atom.name"] == "CA")
        ]

        print(f"The pocket centroid is calculated based on {len(atoms)} CA atoms.")  # TODO logger

        centroid = atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()

        return centroid

    def add_subpocket(
        self, name, anchor_residue_ids, color="blue", anchor_residue_labels=None,
    ):
        """
        Add subpocket based on given anchor residue PDB IDs.

        Parameters
        ----------
        name : str
            Subpocket name.
        anchor_residue_ids : list of (int, str)
            List of anchor residue PDB IDs.
        color : str
            Subpocket color (matplotlib name), blue by default.
        anchor_residue_labels : list of (int, str) or None
            List of anchor residue labels. Must be of same length as anchor_residue_ids.
        """

        subpocket = Subpocket()
        subpocket.from_dataframe(self.data, name, anchor_residue_ids, color, anchor_residue_labels)
        self._subpockets.append(subpocket)

    def add_region(self, name, residue_ids, color="blue", residue_labels=None):
        """
        Add region based on given input residue PDB IDs.

        Parameters
        ----------
        name : str
            Region name.
        residue_ids : list of (int, str)
            List of residue PDB IDs defining the region.
        color : str
            Region color (matplotlib name), blue by default.
        residue_labels : list of (int, str) or None
            List of residue labels. Must be of same length as residue_ids.
        """

        region = Region()
        region.from_dataframe(self.data, name, residue_ids, color, residue_labels)
        self._regions.append(region)

    def visualize(self):
        """
        Visualize the pocket (subpockets, regions, and anchor residues).

        Returns
        -------
        nglview.widget.NGLWidget
            Pocket visualization.
        """

        filepath = str(self.filepath)

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
            color = region["region.color"]
            residue_id = region["residue.id"]
            residue_ngl_ix = residue_id2ix.loc[residue_id]
            scheme_regions_list.append([color, residue_ngl_ix])

        scheme_regions = nglview.color._ColorScheme(scheme_regions_list, label="scheme_regions")
        view.add_representation("cartoon", selection="protein", color=scheme_regions)

        # Show pocket centroids
        view.shape.add_sphere(list(self.centroid), [0, 0, 1], 2, "centroid")

        # Show subpockets
        for index, subpocket in self.subpockets.iterrows():
            center = list(subpocket["subpocket.center"])
            name = subpocket["subpocket.name"]
            color_rgb = colors.to_rgb(subpocket["subpocket.color"])
            view.shape.add_sphere(center, color_rgb, 2, name)

        # Show anchor points
        for index, anchor_residue in self.anchor_residues.iterrows():
            center = list(anchor_residue["anchor_residue.center"])
            color_rgb = colors.to_rgb(anchor_residue["subpocket.color"])
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
        residue_id2ix = self.data[["residue.name", "residue.id"]].drop_duplicates()

        if file_format == "mol2":

            # Map residue names to nglview index (starting from 1)
            residue_id2ix["residue.ngl_ix"] = [str(i) for i in range(1, len(residue_id2ix) + 1)]
            # Cast to Series (PDB IDs as index, NGL index as values)
            residue_id2ix.set_index("residue.id", inplace=True)
            residue_id2ix = residue_id2ix["residue.ngl_ix"]

        else:

            # In this case, PDB ID and nglview index are the same
            residue_id2ix["residue.ngl_ix"] = residue_id2ix["residue.id"]
            residue_id2ix.index = residue_id2ix["residue.id"]
            residue_id2ix = residue_id2ix["residue.ngl_ix"]

        return residue_id2ix
