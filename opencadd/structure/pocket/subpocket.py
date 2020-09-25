"""
opencadd.structure.pocket.subpocket

Defines subpockets.
"""

import numpy as np
import pandas as pd

from .core import Base


class Subpocket(Base):
    """
    Class defining a subpocket.

    Attributes
    ----------
    name : str
        Subpocket name.
    color : str
        Region color name (matplotlib name).
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
            "subpocket.color": [residue.color for residue in self._anchor_residues],
            "anchor_residue.id": [residue.id for residue in self._anchor_residues],
            "anchor_residue.id_alternative": [
                residue.id_alternative for residue in self._anchor_residues
            ],
            "anchor_residue.label": [residue.label for residue in self._anchor_residues],
            "anchor_residue.center": [residue.center for residue in self._anchor_residues],
        }
        anchor_residues_df = pd.DataFrame(anchor_residues_dict)
        anchor_residues_df.insert(0, "subpocket.name", self.name)

        return anchor_residues_df

    def from_dataframe(self, dataframe, name, color, anchor_residue_ids, anchor_residue_labels):
        """
        Set subpocket properties.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural protein data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z".
        name : str
            Subpocket name.
        color : str
            Subpocket color (matplotlib name).
        anchor_residue_ids : list of (int, str)
            List of anchor residue PDB IDs.
        anchor_residue_labels : list of (int, str) or None
            List of anchor residue labels. Must be of same length as residue_ids.
        """

        # Format residue PDB IDs and labels
        anchor_residue_ids, anchor_residue_labels = self._format_residue_ids_and_labels(
            anchor_residue_ids, anchor_residue_labels
        )

        anchor_residues = []

        for residue_id, residue_label in zip(anchor_residue_ids, anchor_residue_labels):

            residue = AnchorResidue()
            residue.from_dataframe(dataframe, residue_id, residue_label, color)
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
            Subpocket color (matplotlib name).
        anchor_residues : list of Residue
            List of anchor residues.
        """

        self.name = name
        self._anchor_residues = anchor_residues
        self.color = color

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


class AnchorResidue(Base):
    """
    Class defining an anchor residue.

    Attributes
    ----------
    id : str
        Residue PDB ID.
    id_alternative : list of str
        Alternative residue PDB ID(s) in case input PDB ID not available.
    label : str
        Residue label, e.g. some non-PDB ID.
    color : str
        Residue color (matplotlib name).
    center : numpy.array
        Coordinates (x, y, z) of the residue center.
    """

    def __init__(self):

        self.id = None
        self.id_alternative = None
        self.label = None
        self.color = None
        self.center = None

    def from_dataframe(self, dataframe, residue_id, residue_label=None, color="green"):
        """
        Set residue properties.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"
        residue_id : str
            Residue PDB ID.
        residue_label : str
            Residue label, e.g. some non-PDB ID (default None).
        color : str
            Residue color (matplotlib name: default green).
        """

        # Set class attributes
        if isinstance(residue_id, int):
            residue_id = str(residue_id)

        self.id = residue_id
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
                self.id_alternative = atoms["residue.id"].to_list()
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
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"

        Returns
        -------
        pandas.DataFrame
            Atom data if DataFrame length is 1, None if length is 0.

        Raises
        ------
        ValueError
            If returned number of atoms is larger than 1.
        """

        atom = dataframe[(dataframe["residue.id"] == self.id) & (dataframe["atom.name"] == "CA")]

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
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"

        Returns
        -------
        pandas.DataFrame
            Atoms data if DataFrame length is 1 or 2, None if length is 0.

        Raises
        ------
        ValueError
            If returned number of atoms is larger than 2.
        """

        residue_id_before = str(int(self.id) - 1)
        residue_id_after = str(int(self.id) + 1)

        atoms = dataframe[
            (dataframe["residue.id"].isin([residue_id_before, residue_id_after]))
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
