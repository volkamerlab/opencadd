"""
opencadd.structure.pocket.api

Defines the opencadd.structure.pocket API.
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

    @classmethod
    def from_file(cls, filepath, name, residue_pdb_ids, residue_labels):
        """
        Initialize Pocket object from structure protein file.

        Attributes
        ----------
        filepath : str or pathlib.Path
            File path to structural protein data.
        name : str
            Name of protein.
        residue_pdb_ids : list of str
            Pocket residue PDB IDs.
        residue_labels : list of str
            Pocket residue labels.

        Returns
        -------
        Pocket
            Pocket object.
        """

        dataframe = DataFrame.from_file(filepath)

        return cls(dataframe, name, residue_pdb_ids, residue_labels)

    @classmethod
    def from_text(cls, text, format, name, residue_pdb_ids, residue_labels):
        """
        Initialize Pocket object from structure file.

        Attributes
        ----------
        text : str
            Structural protein data as string.
        format : str
            Structural protein data format.
        name : str
            Name of protein.
        residue_pdb_ids : list of str
            Pocket residue PDB IDs.
        residue_labels : list of str
            Pocket residue labels.

        Returns
        -------
        Pocket
            Pocket object.
        """

        dataframe = DataFrame.from_text(text, format)

        return cls(dataframe, name, residue_pdb_ids, residue_labels)

    @classmethod
    def from_dataframe(cls, dataframe, name, residue_pdb_ids, residue_labels):
        """
        Initialize Pocket object from DataFrame.

        Attributes
        ----------
        dataframe : pandas.DataFrame
            Structural protein data with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z".
        name : str
            Name of protein.
        residue_pdb_ids : list of str
            Pocket residue PDB IDs.
        residue_labels : list of str
            Pocket residue labels.

        Returns
        -------
        Pocket
            Pocket object.
        """

        return cls(dataframe, name, residue_pdb_ids, residue_labels)

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
                "subpocket.color": [subpocket.color for subpocket in self._subpockets],
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
                    "region.color": [region.color] * n_residues,
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
        self,
        name,
        color,
        anchor_residue_pdb_ids,
        anchor_residue_labels=None,
    ):
        """
        Add subpocket based on given anchor residue PDB IDs.

        Parameters
        ----------
        name : str
            Subpocket name.
        color : str
            Subpocket color (matplotlib name).
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
            Region color (matplotlib name).
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
            color = region["region.color"]
            residue_pdb_id = region["residue.pdb_id"]
            residue_ngl_ix = residue_id2ix.loc[residue_pdb_id]
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
        residue_id2ix = self.data[["residue.name", "residue.pdb_id"]].drop_duplicates()

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
