"""
opencadd.structure.pocket.subpocket

Defines subpockets.
"""

import logging

import numpy as np
import pandas as pd

from .utils import _format_residue_ids_and_labels

_logger = logging.getLogger(__name__)


class Subpocket:
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
        - Anchor residue IDs (user-defined input IDs or alternative
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

    @classmethod
    def from_residue_ids(dataframe, residue_id_ix_mapping, residue_ids, name, color="blue"):
        pass

    @classmethod
    def from_residue_ixs(dataframe, residue_id_ix_mapping, residue_ixs, name, color="blue"):
        pass

    @classmethod
    def from_dataframe(
        cls, dataframe, name, anchor_residue_ids, color="blue", anchor_residue_labels=None
    ):
        """
        Initialize a Subpocket object from a DataFrame.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural protein data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z".
        name : str
            Subpocket name.
        anchor_residue_ids : list of (int, str)
            List of anchor residue IDs.
        color : str
            Subpocket color (matplotlib name), blue by default.
        anchor_residue_labels : list of (int, str) or None
            List of anchor residue labels. Must be of same length as residue_ids.

        Returns
        -------
        opencadd.structure.pocket.subpocket.Subpocket
            Subpocket object.
        """

        # Format residue IDs and labels
        anchor_residue_ids, anchor_residue_labels = _format_residue_ids_and_labels(
            anchor_residue_ids, anchor_residue_labels
        )

        # Get list of AnchorResidue objects
        anchor_residues = []
        for residue_id, residue_label in zip(anchor_residue_ids, anchor_residue_labels):
            residue = AnchorResidue.from_dataframe(dataframe, residue_id, color, residue_label)
            anchor_residues.append(residue)

        subpocket = cls.from_anchor_residues(anchor_residues, name, color)
        return subpocket

    @classmethod
    def from_anchor_residues(cls, anchor_residues, name, color="blue"):
        """
        Initialize a Subpocket object from a list of anchor residues.

        Parameters
        ----------
        anchor_residues : list of Residue
            List of anchor residues.
        name : str
            Subpocket name.
        color : str
            Subpocket color (matplotlib name), blue by default.

        Returns
        -------
        opencadd.structure.pocket.subpocket.Subpocket
            Subpocket object.
        """

        subpocket = cls()
        subpocket.name = name
        subpocket._anchor_residues = anchor_residues
        subpocket.color = color
        subpocket.center = subpocket._centroid()

        return subpocket

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


class AnchorResidue:
    """
    Class defining an anchor residue.

    Attributes
    ----------
    id : str
        Residue ID.
    id_alternative : list of str
        Alternative residue ID(s) in case input ID not available.
    label : str
        Residue label.
    color : str
        Residue color (matplotlib name).
    center : numpy.array
        Coordinates (x, y, z) of the residue center.
    """

    def __init__(self):

        self.id = None
        self.id_alternative = None
        self.label = None
        self.label_alternative = None
        self.color = None
        self.center = None

    @classmethod
    def from_residue_ix(cls, dataframe, residue_id_ix_mapping, residue_ix, color="blue"):
        residue_id = residue_id_ix_mapping.set_index("residue.ix").squeeze().loc[residue_ix]

    @classmethod
    def from_residue_id(cls, dataframe, residue_id_ix_mapping, residue_id, color="blue"):
        """
        Set residue properties.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"
        residue_id : str
            Residue ID.
        color : str
            Residue color (matplotlib name), blue by default.
        residue_label : str
            Residue label (default None).
        """

        residue_ix = residue_id_ix_mapping.set_index("residue.id").squeeze().loc[residue_id]

        # Select atom from residue ID and atom name
        anchor_residue = cls()
        center = anchor_residue._ca_atom_center(dataframe, residue_id)

        # If residue ID exists, get atom coordinates.
        if center:
            return center

        # If not, get atom coordinates for residue before/after.
        else:
            residue_id_before = str(int(self.id) - 1)
            residue_id_after = str(int(self.id) + 1)
            center = anchor_residue._ca_atom_before_after_center(
                dataframe, residue_id_before, residue_id_after
            )
            return center

        if anchor_residue.id_alternative:
            _logger.info(
                f"Anchor residue {anchor_residue.id} is missing, "
                f"use {anchor_residue.id_alternative} instead."
            )
        if anchor_residue.center is None:
            _logger.info(
                f"Anchor residue {anchor_residue.id} and its neighbors are missing "
                f"- no anchor residue could be set."
            )

        return anchor_residue

    def _ca_atom_center(self, dataframe, residue_id):
        """TODO"""

        atom = self._ca_atom(dataframe, residue_id)
        # If residue ID exists, get atom coordinates.
        if atom is not None:
            return atom[["atom.x", "atom.y", "atom.z"]].squeeze().to_numpy()
        else:
            return None

    def _ca_atom_before_after_center(self, dataframe, residue_id_before, residue_id_after):
        """TODO"""

        atoms = anchor_residue._ca_atom_before_after(
            dataframe, residue_id_before, residue_id_after
        )
        if atoms is not None:
            return atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()
        else:
            return None

    def _ca_atom(self, dataframe, residue_id):
        """
        Select a CA atom based on a residue PBD ID.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"
        residue_id : TODO
            TODO

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
            (dataframe["residue.id"] == residue_id) & (dataframe["atom.name"] == "CA")
        ]

        if len(atom) == 1:
            return atom
        elif len(atom) == 0:
            return None
        else:
            raise ValueError(
                f"Unambiguous atom selection. {len(atom)} atoms found instead of 0 or 1."
            )

    def _ca_atom_before_after(self, dataframe, residue_id_before, residue_id_after):
        """
        Select CA atoms from residues before and after a given residue PBD ID.

        Parameters
        ----------
        dataframe : pandas.DataFrame
            Structural data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z"
        residue_id_before : TODO
            TODO
        residue_id_after : TODO
            TODO

        Returns
        -------
        pandas.DataFrame
            Atoms data if DataFrame length is 1 or 2, None if length is 0.

        Raises
        ------
        ValueError
            If returned number of atoms is larger than 2.
        """

        atoms = dataframe[
            (dataframe["residue.id"].isin([residue_id_before, residue_id_after]))
            & (dataframe["atom.name"] == "CA")
        ]

        if len(atoms) == 0:
            return None
        elif len(atoms) <= 2:
            return atoms
        else:
            raise ValueError(
                f"Unambiguous atom selection. {len(atom)} atoms found instead of 0, 1, or 2."
            )
