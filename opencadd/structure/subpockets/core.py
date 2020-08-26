"""
core.py

Defines core classes and functions.
"""

import numpy as np
import pandas as pd


class Pocket:
    """
    Class to define a pocket with 
    - structural data, 
    - subpockets (to be shown as spheres) and 
    - regions (to be highlighted).

    Attributes
    ----------
    TODO
    """

    def __init__(self):

        self.molecule = None
        self._subpockets = None
        self.regions = None

    @property
    def subpockets(self):
        return self._subpockets

    def set_subpockets(
        self, dataframe, names, colors, anchor_residue_ids, anchor_residue_labels,
    ):
        pass


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


class Residue:
    """
    Class to define a residue.

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
        Set residue properties for a given residue PDB ID.

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
        atom = self._atom(dataframe)

        # If residue PDB ID exists, get atom coordinates.
        if atom is not None:
            self.center = atom[["atom.x", "atom.y", "atom.z"]].squeeze().to_numpy()

        # If not, get atom coordinates for residue before/after.
        else:
            atoms = self._atom_before_after(dataframe)
            if atoms is not None:
                self.pdb_id_alternative = atoms["residue.pdb_id"].to_list()
                self.center = atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()
            else:
                self.center = None

    def _atom(self, dataframe):
        """
        Select an atom based on a residue PBD ID and an atom name.

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

    def _atom_before_after(self, dataframe):
        """
        Select atoms with a given atom name for residues before and after a given residue PBD ID.

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
