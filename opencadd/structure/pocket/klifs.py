"""
opencadd.structure.pocket.klifs.KlifsPocket

Defines a KLIFS (kinase) pocket.
"""

from opencadd.databases.klifs import setup_remote
from .core import Pocket


class KlifsPocket(Pocket):
    """
    Extends Pocket to initialize a kinase pocket from a KLIFS structure ID and define standard
    KLIFS regions.
    Refer to Pocket documentation for more information:
    opencadd.structure.pocket.Pocket
    """

    @classmethod
    def from_structure_id(cls, structure_id, subpockets=None):
        """
        Get a KLIFS pocket (remotely by a KLIFS structure ID) that defines the KLIFS regions and
        subpockets.

        Parameters
        ----------
        structure_id : int
            KLIFS structure ID.
        subpockets : pandas.DataFrame
            Subpockets (row) with the following details (columns):
            "anchor_residue.klifs_id" : list of int
                List of anchor residues (KLIFS residue IDs) whose centroid defines the subpocket
                center.
            "subpocket.name" : str
                Subpocket name.
            "subpocket.color" : str
                Subpocket color.

        Returns
        -------
        opencadd.structure.pocket.KlifsPocket
            KLIFS pocket object.
        """

        # Set up remote KLIFS session
        remote = setup_remote()

        # Get pocket and coordinates for a structure (by a KLIFS structure ID)
        pocket = remote.pockets.by_structure_id(structure_id)
        filepath = remote.coordinates.to_pdb(structure_id, ".", entity="complex")

        pocket_3d = cls.from_file(
            filepath,
            pocket["residue.id"].to_list(),
            "example kinase",
            pocket["residue.klifs_id"].to_list(),
        )

        # Add regions
        for (region, color), group in pocket.groupby(
            ["residue.klifs_region_id", "residue.klifs_color"]
        ):
            pocket_3d.add_region(
                region,
                group["residue.id"].to_list(),
                color,
                group["residue.klifs_region_id"].to_list(),
            )

        # Map residue KLIFS IDs > residue ID
        if subpockets is not None:
            subpockets["anchor_residue.ids"] = subpockets["anchor_residue.klifs_ids"].apply(
                lambda x: pocket[pocket["residue.klifs_id"].isin(x)]["residue.id"].to_list()
            )

            # Add subpockets
            for _, subpocket in subpockets.iterrows():
                pocket_3d.add_subpocket(
                    subpocket["subpocket.name"],
                    subpocket["anchor_residue.ids"],
                    subpocket["subpocket.color"],
                    subpocket["anchor_residue.klifs_ids"],
                )

        return pocket_3d
