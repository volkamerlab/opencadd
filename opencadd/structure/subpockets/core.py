"""
core.py

Defines core classes and functions.
"""

import numpy as np
import pandas as pd


class Pocket:
    """
    Class defining a pocket with 
    - structural data (protein/pocket), 
    - subpockets (to be shown as spheres) and 
    - regions (to be highlighted).

    Attributes
    ----------
    data : pandas.DataFrame
        Structural data (protein/pocket) with the following mandatory columns:
        "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z".
    name : str
        Name of protein/pocket.
    _subpockets : list of Subpocket
        List of user-defined subpockets.
    _region : list of Region
        List of user-defined regions.
    """

    def __init__(self, data, name):

        self.data = data
        self.name = name
        self._subpockets = []
        self._regions = []

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

        return subpockets

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
                    "residue.pdb_ids": region.residue_pdb_ids,
                    "residue.label": region.residue_labels,
                }
            )

            regions.append(region)

        regions = pd.concat(regions)

        return regions

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

        return anchor_residues

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
        anchor_residue_labels : list of (int, str)
            List of anchor residue labels.
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
        residue_labels : list of (int, str)
            List of residue labels.
        """

        region = Region()
        region.from_dataframe(self.data, name, color, residue_pdb_ids, residue_labels)
        self._regions.append(region)


class Subpocket:
    """
    Class defining a subpocket.

    Attributes
    ----------
    name : str
        Subpocket name.
    color : str
        Subpocket color.
    center : np.array
        Coordinates (x, y, z) of the subpocket center, 
        i.e. the centroid of all anchor residues' CA atoms.
    _anchor_residues : list of Residue
        List of anchor residues.
    """

    def __init__(self):

        self.name = None
        self.color = None
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
            "anchor_residue.pdb_id": [residue.pdb_id for residue in self._anchor_residues],
            "anchor_residue.pdb_id_alternative": [
                residue.pdb_id_alternative for residue in self._anchor_residues
            ],
            "anchor_residue.label": [residue.label for residue in self._anchor_residues],
            "anchor_residue.center": [residue.center for residue in self._anchor_residues],
        }
        anchor_residues_df = pd.DataFrame(anchor_residues_dict)
        anchor_residues_df.insert(0, "subpocket.name", self.name)
        anchor_residues_df.insert(1, "subpocket.color", self.color)

        return anchor_residues_df

    def from_dataframe(
        self, dataframe, name, color, anchor_residue_pdb_ids, anchor_residue_labels
    ):
        """
        Set subpocket properties.
        
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data (protein/pocket) with the following mandatory columns:
            "residue.pdb_id", "atom.name", "atom.x", "atom.y", "atom.z".
        name : str
            Subpocket name.
        color : str
            Subpocket color.
        anchor_residue_pdb_ids : list of (int, str)
            List of anchor residue PDB IDs.
        anchor_residue_labels : list of (int, str)
            List of anchor residue labels.
        """

        if not anchor_residue_labels:
            anchor_residue_labels = [None] * len(anchor_residue_pdb_ids)

        if len(anchor_residue_pdb_ids) != len(anchor_residue_labels):
            raise ValueError(f"Number of residue PDB IDs and labels must be of same length.")

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

        # Set class attributes
        self.name = name
        self.color = color
        self._anchor_residues = anchor_residues

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


class Region:
    """
    Class defining a region.

    Attributes
    ----------
    name : str
        Region name.
    color : str
        Region color.
    residue_pdb_ids : list of (int, str)
        List of residue PDB IDs defining the region.
    residue_labels : list of (int, str)
        List of residue labels.
    """

    def __init__(self):

        self.name = None
        self.color = None
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
        residue_labels : list of (int, str)
            List of residue labels.
        """

        self.name = name
        self.color = color

        if not residue_labels:
            residue_labels = [None] * len(residue_pdb_ids)

        if len(residue_pdb_ids) != len(residue_labels):
            raise ValueError(f"Number of residue PDB IDs and labels must be of same length.")

        # Cast residue PDB IDs and labels to strings
        residue_pdb_ids = [str(residue_pdb_id) for residue_pdb_id in residue_pdb_ids]
        residue_labels = [str(residue_label) for residue_label in residue_labels]

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


class AnchorResidue:
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
    color : str
        Residue color.
    center : numpy.array
        Coordinates (x, y, z) of the residue center.
    """

    def __init__(self):

        self.pdb_id = None
        self.pdb_id_alternative = None
        self.label = None
        self.color = None
        self.center = None

    def from_dataframe(self, dataframe, residue_pdb_id, residue_label=None, color=None):
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
            Residue color (default None).
        """

        # Set class attributes
        if isinstance(residue_pdb_id, int):
            residue_pdb_id = str(residue_pdb_id)

        self.pdb_id = residue_pdb_id
        self.label = residue_label
        self.color = color

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
