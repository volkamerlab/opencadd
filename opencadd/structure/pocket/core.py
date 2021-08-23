"""
opencadd.structure.pocket.core

Defines pockets.
"""

import logging
from pathlib import Path

import pandas as pd
from opencadd.io import DataFrame

from .base import PocketBase
from .features.region import Region
from .features.subpocket import Subpocket
from .features.anchor import AnchorResidue

_logger = logging.getLogger(__name__)


class Pocket(PocketBase):
    """
    Class defining a pocket with
    - structural protein data,
    - subpockets (to be shown as spheres) and
    - regions (to be highlighted).

    Attributes
    ----------
    name : str
        Name of protein.
    data : pd.DataFrame
        Structural protein data with the following mandatory columns:
        "residue.id", "atom.name", "atom.x", "atom.y", "atom.z".
    center : numpy.array
        Pocket center (coordinates).
    _text : str
        Structural protein data as string (file content).
    _extension : str
        Structural protein data format (file extension).
    _residue_ids : list of int
        Pocket residue IDs.
    _residue_ixs : list of (int or None)
        Pocket residue indices.
    _subpockets : list of Subpocket
        List of user-defined subpockets.
    _region : list of Region
        List of user-defined regions.

    Properties
    ----------
    residues
    subpockets
    regions
    anchor_residues
    """

    def __init__(self):

        self.name = None
        self.data = None
        self.center = None
        self._text = None
        self._extension = None
        self._residue_ids = None
        self._residue_ixs = None
        self._subpockets = []
        self._regions = []

    @classmethod
    def from_file(cls, filepath, residue_ids, residue_ixs=None, name=None):
        """
        Initialize Pocket object from structure protein file.

        Attributes
        ----------
        filepath : str or pathlib.Path
            File path to structural protein data.
        residue_ids : list of int
            Pocket residue IDs.
        residue_ixs : None or list of int
            Pocket residue indices. Set to None by default.
        name : str or None
            Name of protein (default: None).

        Returns
        -------
        opencadd.structure.pocket.Pocket
            Pocket object.
        """

        filepath = Path(filepath)
        extension = filepath.suffix[1:]
        with open(filepath, "r") as f:
            text = f.read()
        pocket = cls.from_text(text, extension, residue_ids, residue_ixs, name)
        return pocket

    @classmethod
    def from_text(cls, text, extension, residue_ids, residue_ixs=None, name=None):
        """
        Initialize Pocket object from structure protein text.

        Attributes
        ----------
        text : str
            Structural protein data as string (file content).
        extension : str
            Structural protein data format (file extension).
        residue_ids : list of int
            Pocket residue IDs.
        residue_ixs : None or list of int
            Pocket residue indices. Set to None by default.
        name : str or None
            Name of protein (default: None).

        Returns
        -------
        opencadd.structure.pocket.Pocket
            Pocket object.
        """

        pocket = cls()
        pocket.name = name
        pocket._text = text
        pocket._extension = extension
        pocket._residue_ids, pocket._residue_ixs = pocket._format_residue_ids_and_ixs(
            residue_ids, residue_ixs, "set pocket residues"
        )
        pocket.data = pocket._set_data()
        pocket.center = pocket._set_center()
        return pocket

    def _set_data(self):
        """
        Load atoms as DataFrames from text.
        Keep only atoms with int-castable residue PDB IDs.

        Returns
        -------
        pd.DataFrame
            Structural protein data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z".
        """

        # Load atoms as DataFrame from text
        # Note: Column "residue.id" is of string type
        dataframe = DataFrame.from_text(self._text, self._extension)

        # Cast the IDs to str, so that they can match the DataFrame's ID
        residue_ids = [str(residue_id) for residue_id in self._residue_ids if residue_id]
        # Fetch all atoms matching residues IDs and cast them back to integers
        # This step will remove any residues with insertion codes,
        # e.g. "612A" in structure 7oam (PDB ID)
        dataframe = dataframe[dataframe["residue.id"].isin(residue_ids)]
        dataframe = dataframe.astype({"residue.id": "int32"})

        return dataframe

    def _set_center(self):
        """
        Pocket center, i.e. the centroid of all input residues' CA atoms.

        Returns
        ----------
        numpy.array
            Pocket center (coordinates).
        """

        dataframe = self.data

        atoms = dataframe[
            (dataframe["residue.id"].isin(self._residue_ids)) & (dataframe["atom.name"] == "CA")
        ]

        if len(atoms) != len(self._residue_ids):
            _logger.info(
                f"Pocket {self.name}: Missing pocket CA atoms. "
                f"The pocket center is calculated based on {len(atoms)} CA atoms "
                f"(total number of pocket residues is {len(self._residue_ids)})."
            )

        center = atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()

        return center

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

        subpockets = pd.DataFrame([subpocket.subpocket for subpocket in self._subpockets])
        return subpockets.reset_index(drop=True)

    @property
    def regions(self):
        """
        All pocket's regions.

        Returns
        -------
        pandas.DataFrame
            Name, color, involved residue IDs and indices (columns) for all regions.
        """
        if self._regions == []:
            return None

        regions = [region.region for region in self._regions]
        regions = pd.concat(regions)
        regions = regions.astype({"residue.id": "Int32"})
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
            - Anchor residue IDs (user-defined input IDs or alternative
            IDs if input was not available)
            - Anchor residue indices
            - The anchor residue centers (coordinates)
        """

        if self._subpockets == []:
            return None

        anchor_residues = [subpocket.anchor_residues for subpocket in self._subpockets]
        anchor_residues = pd.concat(anchor_residues)
        anchor_residues = anchor_residues.astype({"anchor_residue.id": "Int32"})

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

    def add_subpocket(
        self,
        name,
        anchor_residue_ids=None,
        anchor_residue_ixs=None,
        color="blue",
    ):
        """
        Add subpocket based on given anchor residue IDs.

        Parameters
        ----------
        name : str
            Subpocket name.
        anchor_residue_ids : list of (int or str) or None
            List of anchor residue IDs.
            Note: If "anchor_residue_ids" and "anchor_residue_ix" are both set,
            "anchor_residue_ix" will be used.
        anchor_residue_ixs : list of (int or str) or None
            List of anchor residue indices.
            Note: If "anchor_residue_ids" and "anchor_residue_ix" are both set,
            "anchor_residue_ix" will be used.
        color : str
            Subpocket color (matplotlib name), blue by default.
        """

        if anchor_residue_ids and anchor_residue_ixs:
            raise ValueError(f"Please set only anchor residue PDB IDs or indices - not both.")
        if anchor_residue_ixs:
            subpocket = self._subpocket_by_residue_ixs(anchor_residue_ixs, name, color)
        else:
            subpocket = self._subpocket_by_residue_ids(anchor_residue_ids, name, color)
        self._subpockets.append(subpocket)

    def add_region(self, name, residue_ids=None, residue_ixs=None, color="blue"):
        """
        Add region based on given input residue IDs.

        Parameters
        ----------
        name : str
            Region name.
        residue_ids : list of (int or str) or None
            List of residue IDs defining the region.
            Note: If "residue_ids" and "residue_ix" are both set, "residue_ix" will be used.
        residue_ixs : list of (int or str) or None
            List of residue indices.
            Note: If "residue_ids" and "residue_ix" are both set, "residue_ix" will be used.
        color : str
            Region color (matplotlib name), blue by default.
        """

        if residue_ids and residue_ixs:
            raise ValueError(f"Please set only residue PDB IDs or indices - not both.")
        if residue_ixs:
            residue_ids = [self._residue_ix2id(residue_ix) for residue_ix in residue_ixs]
        else:
            residue_ixs = [self._residue_id2ix(residue_id) for residue_id in residue_ids]

        residue_ids, residue_ixs = self._format_residue_ids_and_ixs(
            residue_ids, residue_ixs, "set region residues"
        )
        region = Region(name, residue_ids, residue_ixs, color)
        self._regions.append(region)

    @property
    def ca_atoms(self):
        """
        CA atoms of the pocket residues.

        Returns
        -------
        pandas.DataFrame
            Structural data for CA atoms of the pocket residues.
        """
        return self._ca_atoms(*self._residue_ids)

    def _ca_atoms(self, *residue_ids):
        r"""
        Select a CA atoms based on residue PBD ID(s).

        Parameters
        ----------
        \*residue_ids : int
            Residue PDB ID(s).

        Returns
        -------
        pandas.DataFrame or None
            Structural data for CA atoms of input residues.
            If no CA atoms available, returns None.

        Raises
        ------
        ValueError
            If returned number of CA atoms is larger than 1.
        """

        ca_atoms = self.data[  # pylint: disable=E1136
            (self.data["residue.id"].isin(residue_ids))  # pylint: disable=E1136
            & (self.data["atom.name"] == "CA")  # pylint: disable=E1136
        ]

        if len(ca_atoms) <= len(residue_ids):
            return ca_atoms
        else:
            raise ValueError(
                f"More CA atoms ({len(ca_atoms)}) found than input residues ({len(residue_ids)})."
            )

    def _ca_atoms_center(self, *residue_ids):
        r"""
        Calculate centroid of all CA atoms based on residue PBD ID(s).

        Parameters
        ----------
        \*residue_ids : str
            Residue PDB ID(s).

        Returns
        -------
        numpy.array or None
            Centroid of all CA atoms of input residues.
            If no CA atoms available, returns None.
        """

        ca_atoms = self._ca_atoms(*residue_ids)

        # If residue ID exists, get atom coordinates.
        if len(ca_atoms) == 1:
            center = ca_atoms[["atom.x", "atom.y", "atom.z"]].squeeze().to_numpy()
        elif len(ca_atoms) > 1:
            center = ca_atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()
        else:
            center = None

        ca_atoms_residue_ids = ca_atoms["residue.id"].to_list()
        if len(ca_atoms_residue_ids) == 0:
            ca_atoms_residue_ids = None

        return ca_atoms_residue_ids, center

    def _anchor_residue_by_residue_id(self, residue_id, color="blue", subpocket_name=None):
        """
        Get anchor residue (AnchorResidue object) based on a selected residue PDB ID.

        Parameters
        ----------
        residue_id : int
            Residue PDB ID.
        color : str
            Color name (matplotlib).
        subpocket_name : str or None
            Subpocket name.

        Returns
        -------
        opencadd.structure.pocket.AnchorResidue
            Anchor residue.
        """

        residue_id_alternative = None
        _, center = self._ca_atoms_center(residue_id)

        if center is None:
            residue_id_before = residue_id - 1
            residue_id_after = residue_id + 1
            residue_ids = [residue_id_before, residue_id_after]
            residue_id_alternative, center = self._ca_atoms_center(*residue_ids)

        subpocket_anchor = AnchorResidue(
            center, residue_id, residue_id_alternative, None, color, subpocket_name, self.name
        )

        return subpocket_anchor

    def _anchor_residue_by_residue_ix(self, residue_ix, color="blue", subpocket_name=None):
        """
        Get anchor residue (AnchorResidue object) based on a selected residue index.

        Parameters
        ----------
        residue_ix : int or str
            Residue index.
        color : str
            Color name (matplotlib).
        subpocket_name : str or None
            Subpocket name.

        Returns
        -------
        opencadd.structure.pocket.AnchorResidue
            Anchor residue.
        """

        residue_id = self._residue_ix2id(residue_ix)
        residue_id_alternative = None

        if residue_id:
            subpocket_anchor = self._anchor_residue_by_residue_id(
                residue_id, color, subpocket_name
            )
            subpocket_anchor.residue_ix = residue_ix
        else:
            # Get residue indices before and after
            residue_ix_before = residue_ix - 1
            residue_ix_after = residue_ix + 1
            # Get corresponding residue IDs
            residue_id_before = self._residue_ix2id(residue_ix_before)
            residue_id_after = self._residue_ix2id(residue_ix_after)
            residue_ids = [residue_id_before, residue_id_after]
            residue_ids = [residue_id for residue_id in residue_ids if residue_id is not None]
            # Add one more check: Residue IDs can only be separated by one residue!
            if len(residue_ids) == 2:
                if int(residue_ids[1]) - int(residue_ids[0]) != 2:
                    residue_ids = []
            # Get center
            residue_id_alternative, center = self._ca_atoms_center(*residue_ids)

            subpocket_anchor = AnchorResidue(
                center,
                residue_id,
                residue_id_alternative,
                residue_ix,
                color,
                subpocket_name,
                self.name,
            )

        return subpocket_anchor

    def _subpocket_by_residue_ids(self, residue_ids, name=None, color="blue"):
        """
        Generate subpocket based on residue PDB ID(s).

        Parameters
        ----------
        residue_ids : list of (int or str)
            Residue PDB ID(s).
        name : str or None
            Subpocket name.
        color : str
            Color name (matplotlib).

        Returns
        -------
        opencadd.structure.pocket.Subpocket
            Subpocket.
        """

        anchor_residues = []
        for residue_id in residue_ids:
            anchor_residue = self._anchor_residue_by_residue_id(residue_id, color, name)
            anchor_residues.append(anchor_residue)

        subpocket = Subpocket(anchor_residues, name, color)
        return subpocket

    def _subpocket_by_residue_ixs(self, residue_ixs, name=None, color="blue"):
        """
        Generate subpocket based on residue indices.

        Parameters
        ----------
        residue_ixs : list of (int or str)
            Residue indices.
        name : str or None
            Subpocket name.
        color : str
            Color name (matplotlib).

        Returns
        -------
        opencadd.structure.pocket.Subpocket
            Subpocket.
        """

        anchor_residues = []
        for residue_ix in residue_ixs:
            anchor_residue = self._anchor_residue_by_residue_ix(residue_ix, color, name)
            anchor_residues.append(anchor_residue)

        subpocket = Subpocket(anchor_residues, name, color)
        return subpocket
