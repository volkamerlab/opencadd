"""
opencadd.structure.pocket.core

Defines pockets.
"""

import logging
from pathlib import Path

import pandas as pd

from opencadd.io import DataFrame
from .region import Region
from .subpocket import Subpocket
from .anchor import AnchorResidue
from .utils import _format_residue_ids_and_ixs

_logger = logging.getLogger(__name__)


class Pocket:
    """
    Class defining a pocket with
    - structural protein data,
    - subpockets (to be shown as spheres) and
    - regions (to be highlighted).

    Attributes
    ----------
    data
    residues
    subpockets
    regions
    anchor_residues
    center
    name : str
        Name of protein.
    _text : str
        Structural protein data as string (file content).
    _extension : str
        Structural protein data format (file extension).
    _residue_ids : list of str
        Pocket residue IDs.
    _residue_ixs : list of str
        Pocket residue indices.
    _subpockets : list of Subpocket
        List of user-defined subpockets.
    _region : list of Region
        List of user-defined regions.
    """

    def __init__(self):

        self.name = None
        self.data = None
        self._text = None
        self._extension = None
        self._residue_ids = None
        self._residue_ixs = None
        self._subpockets = []
        self._regions = []

    @classmethod
    def from_file(cls, filepath, residue_ids, name="", residue_ixs=None):
        """
        Initialize Pocket object from structure protein file.

        Attributes
        ----------
        filepath : str or pathlib.Path
            File path to structural protein data.
        residue_ids : list of str
            Pocket residue IDs.
        name : str
            Name of protein (default: empty string).
        residue_ixs : None or list of str
            Pocket residue indices. Set to None by default.

        Returns
        -------
        opencadd.structure.pocket.Pocket
            Pocket object.
        """

        filepath = Path(filepath)
        extension = filepath.suffix[1:]
        with open(filepath, "r") as f:
            text = f.read()
        pocket = cls.from_text(text, extension, residue_ids, name, residue_ixs)
        return pocket

    @classmethod
    def from_text(cls, text, extension, residue_ids, name="", residue_ixs=None):
        """
        Initialize Pocket object from structure protein text.

        Attributes
        ----------
        text : str
            Structural protein data as string (file content).
        extension : str
            Structural protein data format (file extension).
        residue_ids : list of str
            Pocket residue IDs.
        name : str
            Name of protein (default: empty string).
        residue_ixs : None or list of str
            Pocket residue indices. Set to None by default.

        Returns
        -------
        opencadd.structure.pocket.Pocket
            Pocket object.
        """

        dataframe = DataFrame.from_text(text, extension)
        pocket = cls._from_dataframe(dataframe, residue_ids, name, residue_ixs)
        pocket._text = text
        pocket._extension = extension
        return pocket

    @classmethod
    def _from_dataframe(cls, dataframe, residue_ids, name="", residue_ixs=None):
        """
        Initialize Pocket object from structure DataFrame.

        Attributes
        ----------
        dataframe : pd.DataFrame
            Structural protein data with the following mandatory columns:
            "residue.id", "atom.name", "atom.x", "atom.y", "atom.z".
        residue_ids : list of str
            Pocket residue IDs.
        name : str
            Name of protein (default: empty string).
        residue_ixs : None or list of str
            Pocket residue indices. Set to None by default.

        Returns
        -------
        opencadd.structure.pocket.Pocket
            Pocket object.
        """

        pocket = cls()
        pocket.name = name
        pocket._residue_ids, pocket._residue_ixs = _format_residue_ids_and_ixs(
            residue_ids, residue_ixs
        )
        pocket.data = dataframe[dataframe["residue.id"].isin(pocket._residue_ids)]
        return pocket

    @property
    def residues(self):
        """
        All pocket's residues.

        Returns
        -------
        pandas.DataFrame
            Residue ID and residue index (columns) for all pocket residues (rows).
        """

        residues = {"residue.id": self._residue_ids, "residue.ix": self._residue_ixs}
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

        regions = []

        for region in self._regions:

            n_residues = len(region.residue_ids)

            region = pd.DataFrame(
                {
                    "region.name": [region.name] * n_residues,
                    "region.color": [region.color] * n_residues,
                    "residue.id": region.residue_ids,
                    "residue.ix": region.residue_ixs,
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
            - Anchor residue IDs (user-defined input IDs or alternative
            IDs if input was not available)
            - Anchor residue indices
            - The anchor residue centers (coordinates)
        """

        if self._subpockets == []:
            return None

        anchor_residues = [subpocket.anchor_residues for subpocket in self._subpockets]
        anchor_residues = pd.concat(anchor_residues)

        return anchor_residues.reset_index(drop=True)

    @property
    def center(self):
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
                f"Missing pocket CA atoms. "
                f"The pocket center is calculated based on {len(atoms)} CA atoms "
                f"(total number of pocket residues is {len(self._residue_ids)})."
            )

        center = atoms[["atom.x", "atom.y", "atom.z"]].mean().to_numpy()

        return center

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
        anchor_residue_ids,
        color="blue",
        anchor_residue_ixs=None,
    ):
        """
        Add subpocket based on given anchor residue IDs.

        Parameters
        ----------
        name : str
            Subpocket name.
        anchor_residue_ids : list of (int, str)
            List of anchor residue IDs.
        color : str
            Subpocket color (matplotlib name), blue by default.
        anchor_residue_ixs : list of (int, str) or None
            List of anchor residue indices. Must be of same length as anchor_residue_ids.
        """

        if anchor_residue_ixs:
            subpocket = self._subpocket_by_residue_ixs(anchor_residue_ixs, name, color)
        else:
            subpocket = self._subpocket_by_residue_ids(anchor_residue_ids, name, color)
        self._subpockets.append(subpocket)

    def add_region(self, name, residue_ids, color="blue", residue_ixs=None):
        """
        Add region based on given input residue IDs.

        Parameters
        ----------
        name : str
            Region name.
        residue_ids : list of (int, str)
            List of residue IDs defining the region.
        color : str
            Region color (matplotlib name), blue by default.
        residue_ixs : list of (int, str) or None
            List of residue indices. Must be of same length as residue_ids.
        """

        region = Region()
        region.from_dataframe(self.data, name, residue_ids, color, residue_ixs)
        self._regions.append(region)

    def _residue_ix2id(self, residue_ix):
        """
        Get residue PDB ID from residue index.

        Parameters
        ----------
        residue_ix : int or str
            Residue index.

        Returns
        -------
        str
            Residue PDB ID.
        """

        residue_ix = str(residue_ix)
        residues = self.residues
        residues = residues[~residues["residue.ix"].isin(["_", "-", "", " ", None])]
        try:
            residue_id = self.residues.set_index("residue.ix").squeeze().loc[residue_ix]
        except KeyError:
            residue_id = None
        return residue_id

    def _residue_id2ix(self, residue_id):
        """
        Get residue index from residue PDB ID.

        Parameters
        ----------
        residue_id : int or str
            Residue PDB ID.

        Returns
        -------
        str
            Residue index.
        """

        residue_id = str(residue_id)
        residues = self.residues
        residues = residues[~residues["residue.id"].isin(["_", "-", "", " ", None])]
        try:
            residue_ix = residues.set_index("residue.id").squeeze().loc[residue_id]
        except KeyError:
            residue_ix = None
        return residue_ix

    def _ca_atoms(self, *residue_ids):
        r"""
        Select a CA atoms based on residue PBD ID(s).

        Parameters
        ----------
        \*residue_ids : str
            Residue PDB ID(s).

        Returns
        -------
        pandas.DataFrame
            Atom data if DataFrame length is 1, None if length is 0. TODO

        Raises
        ------
        ValueError
            If returned number of atoms is larger than 1. TODO
        """

        residue_ids = [str(residue_id) for residue_id in residue_ids]
        ca_atoms = self.data[
            (self.data["residue.id"].isin(residue_ids)) & (self.data["atom.name"] == "CA")
        ]

        if len(ca_atoms) <= len(residue_ids):
            return ca_atoms
        else:
            raise ValueError(
                f"More CA atoms ({len(ca_atoms)}) found than input residues ({len(residue_ids)})."
            )

    def _ca_atoms_center(self, *residue_ids):
        """TODO"""

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

    def _anchor_residue_by_residue_id(self, residue_id, color="blue"):
        """
        Get anchor residue (AnchorResidue object) based on a selected residue PDB ID.

        Parameters
        ----------
        residue_id : int or str
            Residue PDB ID.
        color : str
            Color name (matplotlib).

        Returns
        -------
        opencadd.structure.pocket.AnchorResidue
            Anchor residue.
        """

        residue_id = str(residue_id)
        residue_id_alternative = None
        _, center = self._ca_atoms_center(residue_id)

        if center is None:
            residue_id_before = str(int(residue_id) - 1)
            residue_id_after = str(int(residue_id) + 1)
            residue_ids = [residue_id_before, residue_id_after]
            residue_id_alternative, center = self._ca_atoms_center(*residue_ids)

        subpocket_anchor = AnchorResidue(center, residue_id, residue_id_alternative, None, color)

        return subpocket_anchor

    def _anchor_residue_by_residue_ix(self, residue_ix, color="blue"):
        """
        Get anchor residue (AnchorResidue object) based on a selected residue index.

        Parameters
        ----------
        residue_ix : int or str
            Residue index.
        color : str
            Color name (matplotlib).

        Returns
        -------
        opencadd.structure.pocket.AnchorResidue
            Anchor residue.
        """

        residue_ix = str(residue_ix)
        residue_id = self._residue_ix2id(residue_ix)
        residue_id_alternative = None

        if residue_id:
            subpocket_anchor = self._anchor_residue_by_residue_id(residue_id, color)
            subpocket_anchor.residue_ix = residue_ix
        else:
            # Get residue indices before and after
            residue_ix_before = str(int(residue_ix) - 1)
            residue_ix_after = str(int(residue_ix) + 1)
            # Get corresponding residue IDs
            residue_id_before = self._residue_ix2id(residue_ix_before)
            residue_id_after = self._residue_ix2id(residue_ix_after)
            residue_ids = [residue_id_before, residue_id_after]
            residue_ids = [residue_id for residue_id in residue_ids if residue_id is not None]
            # Get center
            residue_id_alternative, center = self._ca_atoms_center(*residue_ids)

            subpocket_anchor = AnchorResidue(
                center, residue_id, residue_id_alternative, residue_ix, color
            )

        return subpocket_anchor

    def _subpocket_by_residue_ids(self, residue_ids, name=None, color="blue"):
        """TODO"""

        anchor_residues = []
        for residue_id in residue_ids:
            anchor_residue = self._anchor_residue_by_residue_id(residue_id, color)
            anchor_residues.append(anchor_residue)

        subpocket = Subpocket(anchor_residues, name, color)
        return subpocket

    def _subpocket_by_residue_ixs(self, residue_ixs, name=None, color="blue"):
        """TODO"""

        anchor_residues = []
        for residue_ix in residue_ixs:
            anchor_residue = self._anchor_residue_by_residue_ix(residue_ix, color)
            anchor_residues.append(anchor_residue)

        subpocket = Subpocket(anchor_residues, name, color)
        return subpocket
