"""
Tests for opencadd.structure.pocket.core
"""

from pathlib import Path

import pandas as pd
import pytest

from opencadd.structure.pocket import Pocket, AnchorResidue

PATH_TEST_DATA = Path(__name__).parent / "opencadd/tests/data/pocket"


class TestsPocket:
    """
    Test Pocket class methods.
    """

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, name, data_residue_ids",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 13],
                None,
                "example kinase",
                [12, 13]
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                ["15", "16", "_"],
                [4, 5, 6],
                "example kinase",
                [15, 16]
            ),
        ],
    )
    def test_from_file(self, filepath, residue_ids, residue_ixs, name, data_residue_ids):
        """
        Initialize class from file and test class attributes and properties.
        """

        pocket = Pocket.from_file(filepath, residue_ids, residue_ixs, name)

        # Test attributes
        assert pocket.name == name
        assert isinstance(pocket._text, str)
        assert pocket._extension == filepath.suffix[1:]
        assert isinstance(pocket.data, pd.DataFrame)
        assert list(pocket.data["residue.id"].unique()) == data_residue_ids
        # Note: The following attributes/properties are tested for parent class BasePocket:
        # - _residue_ids
        # - _residue_ixs
        # - residues

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, pocket_residue_ixs, region_name, region_color",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                [16, 47, 80],
                "hinge",
                "magenta",
            )
        ],
    )
    def test_add_regions(
        self,
        filepath,
        pocket_residue_ids,
        pocket_residue_ixs,
        region_name,
        region_color,
    ):
        """
        Test adding regions and associated regions property and attribute.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids, pocket_residue_ixs)
        assert pocket.regions == None

        # Add region
        pocket.add_region(region_name, residue_ids=pocket_residue_ids, color=region_color)
        pocket.add_region(region_name, residue_ixs=pocket_residue_ixs, color=region_color)

        # Test the subpocket in the _subpockets attribute
        assert pocket._regions[0].residue_ids == pocket._regions[1].residue_ids
        assert pocket._regions[0].residue_ixs == pocket._regions[1].residue_ixs

        # Test regions property
        assert isinstance(pocket.regions, pd.DataFrame)
        assert pocket.regions.columns.to_list() == [
            "region.name",
            "region.color",
            "residue.id",
            "residue.ix",
        ]

        # Test an example region in the _regions attribute
        region = pocket._regions[0]
        assert region.name == region_name
        assert region.color == region_color

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, pocket_residue_ixs, subpocket_name, subpocket_color, subpocket_center",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [73, 128, 193],
                [16, 47, 80],
                "hinge_region",
                "magenta",
                [1.957633, 21.923767, 41.69003333],
            )
        ],
    )
    def test_add_subpocket(
        self,
        filepath,
        pocket_residue_ids,
        pocket_residue_ixs,
        subpocket_name,
        subpocket_color,
        subpocket_center,
    ):
        """
        Test adding subpocket and associated subpockets property and attribute.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids, pocket_residue_ixs)
        assert pocket.subpockets == None
        assert pocket.anchor_residues == None

        # Add subpocket
        pocket.add_subpocket(
            subpocket_name, anchor_residue_ids=pocket_residue_ids, color=subpocket_color
        )
        pocket.add_subpocket(
            subpocket_name, anchor_residue_ixs=pocket_residue_ixs, color=subpocket_color
        )

        # Test the subpocket in the _subpockets attribute
        assert all(pocket._subpockets[0].center == pocket._subpockets[1].center)
        subpocket = pocket._subpockets[0]
        assert subpocket.name == subpocket_name
        assert subpocket.color == subpocket_color
        assert pytest.approx(subpocket.center, abs=1.0e-6) == subpocket_center

        # Test subpockets property
        assert isinstance(pocket.subpockets, pd.DataFrame)
        assert pocket.subpockets.columns.to_list() == [
            "subpocket.name",
            "subpocket.color",
            "subpocket.center",
        ]

        # Test anchor_residues property
        assert isinstance(pocket.anchor_residues, pd.DataFrame)
        assert pocket.anchor_residues.columns.to_list() == [
            "subpocket.name",
            "anchor_residue.color",
            "anchor_residue.id",
            "anchor_residue.id_alternative",
            "anchor_residue.ix",
            "anchor_residue.center",
        ]

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, subpocket_residue_ids, subpocket_center, subpocket_name, subpocket_color",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [16, 51, 53],
                [16, 51],
                [11.512, 17.917, 33.356],
                "subpocket1",
                "red",
            ),
        ],
    )
    def test_subpocket_by_residue_ids(
        self,
        filepath,
        pocket_residue_ids,
        subpocket_residue_ids,
        subpocket_center,
        subpocket_name,
        subpocket_color,
    ):
        """
        Test subpocket definitions generated based on residue PDB IDs.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids)
        subpocket = pocket._subpocket_by_residue_ids(
            subpocket_residue_ids, subpocket_name, subpocket_color
        )

        assert subpocket.name == subpocket_name
        assert subpocket.color == subpocket_color
        for i, j in zip(subpocket.center, subpocket_center):
            assert pytest.approx(i, abs=1.0e-3) == j

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, pocket_residue_ixs, subpocket_residue_ixs, subpocket_center, subpocket_name, subpocket_color",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [16, 51, 53],
                [5, 20, 22],
                [5, 20],
                [11.512, 17.917, 33.356],
                "subpocket1",
                "red",
            ),
        ],
    )
    def test_subpocket_by_residue_ixs(
        self,
        filepath,
        pocket_residue_ids,
        pocket_residue_ixs,
        subpocket_residue_ixs,
        subpocket_center,
        subpocket_name,
        subpocket_color,
    ):
        """
        Test subpocket definitions generated based on residue indices.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids, pocket_residue_ixs)
        subpocket = pocket._subpocket_by_residue_ixs(
            subpocket_residue_ixs, subpocket_name, subpocket_color
        )

        assert subpocket.name == subpocket_name
        assert subpocket.color == subpocket_color
        for i, j in zip(subpocket.center, subpocket_center):
            assert pytest.approx(i, abs=1.0e-3) == j

    @pytest.mark.parametrize(
        "filepath, residue_ids, anchor_residue_id, anchor_residue_id_alternative, anchor_residue_color, anchor_residue_center",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [12, 16, 85, 87],
                12,
                None,
                "yellow",
                [8.1173, 16.3077, 51.996],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                86,
                [85, 87],
                "yellow",
                [3.8648999, 26.4354, 42.65935],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                17,
                [16],
                "yellow",
                [9.3363, 11.0014, 42.1329],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                18,
                None,
                "yellow",
                None,
            ),
        ],
    )
    def test_anchor_residue_by_residue_id(
        self,
        filepath,
        residue_ids,
        anchor_residue_id,
        anchor_residue_id_alternative,
        anchor_residue_color,
        anchor_residue_center,
    ):
        """
        Test anchor residue definitions generated based on residue PDB IDs.
        """

        pocket = Pocket.from_file(filepath, residue_ids)
        anchor_residue = pocket._anchor_residue_by_residue_id(
            anchor_residue_id, anchor_residue_color
        )
        assert isinstance(anchor_residue, AnchorResidue)
        assert anchor_residue.residue_id == anchor_residue_id
        assert anchor_residue.residue_id_alternative == anchor_residue_id_alternative
        assert anchor_residue.residue_ix == None
        assert anchor_residue.color == anchor_residue_color
        if anchor_residue_center:
            assert pytest.approx(anchor_residue.center[0], abs=1.0e-3) == anchor_residue_center[0]

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, anchor_residue_ix, anchor_residue_id, anchor_residue_id_alternative, anchor_residue_color, anchor_residue_center",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                1,
                12,
                None,
                "yellow",
                [8.1173, 16.3077, 51.996],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                45,
                None,
                [85, 87],
                "yellow",
                [3.8648999, 26.4354, 42.65935],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                6,
                None,
                [16],
                "yellow",
                [9.3363, 11.0014, 42.1329],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                7,
                None,
                None,
                "yellow",
                None,
            ),
        ],
    )
    def test_anchor_residue_by_residue_ix(
        self,
        filepath,
        residue_ids,
        residue_ixs,
        anchor_residue_ix,
        anchor_residue_id,
        anchor_residue_id_alternative,
        anchor_residue_color,
        anchor_residue_center,
    ):
        """
        Test anchor residue definitions generated based on residue indices.
        """

        pocket = Pocket.from_file(filepath, residue_ids, residue_ixs)
        anchor_residue = pocket._anchor_residue_by_residue_ix(
            anchor_residue_ix, anchor_residue_color
        )
        assert isinstance(anchor_residue, AnchorResidue)
        assert anchor_residue.residue_ix == anchor_residue_ix
        if anchor_residue_id:
            assert anchor_residue.residue_id == anchor_residue_id
        assert anchor_residue.residue_id_alternative == anchor_residue_id_alternative
        assert anchor_residue.color == anchor_residue_color
        if anchor_residue_center:
            assert pytest.approx(anchor_residue.center[0], abs=1.0e-3) == anchor_residue_center[0]

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, residue_ids, n_ca_atoms, center",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [12, 16, 85, 87],
                [12],
                1,
                [8.117, 16.308, 51.996],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [12, 16],
                2,
                [8.727, 13.655, 47.064],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [13],
                0,
                None,
            ),
        ],
    )
    def test_ca_atoms_and_ca_atoms_center(
        self, filepath, pocket_residue_ids, residue_ids, n_ca_atoms, center
    ):
        """
        Test CA atoms retrieval and center calculation based on a set of residue PDB IDs.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids)

        # Test method _ca_atoms()
        ca_atoms = pocket._ca_atoms(*residue_ids)
        assert len(ca_atoms) == n_ca_atoms

        # Test method _ca_atoms_center()
        _, ca_atoms_center = pocket._ca_atoms_center(*residue_ids)
        if center:
            for i, j in zip(ca_atoms_center, center):
                assert pytest.approx(i, abs=1.0e-3) == j

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, ca_atoms_residue_ids, pocket_center",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",  # KLIFS ID 3834
                [127, 128, 129],
                [127, 128, 129],
                [-1.178, 23.860, 45.092],
            ),
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",  # KLIFS ID 3834
                [1, 127, 128, 129],
                [127, 128, 129],
                [-1.178, 23.860, 45.092],
            ),
        ],
    )
    def test_center_and_ca_atoms(self, filepath, pocket_residue_ids, ca_atoms_residue_ids, pocket_center):
        """
        Test the pocket CA atoms and center calculation based on the pocket's residue PDB IDs.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids)


        # Test property ca_atom
        assert len(pocket.ca_atoms) == len(ca_atoms_residue_ids)

        # Test property center
        for i, j in zip(pocket.center, pocket_center):
            assert pytest.approx(i, abs=1.0e-3) == j

    @pytest.mark.parametrize(
        "residue_ids, residue_ixs, residue_ids_formatted, residue_ixs_formatted",
        [
            ([1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3]),
            ([1, 2], None, [1, 2], [None, None]),
        ],
    )
    def test_format_residue_ids_and_ixs(
        self, residue_ids, residue_ixs, residue_ids_formatted, residue_ixs_formatted
    ):

        pocket = Pocket()
        residue_ids, residue_ixs = pocket._format_residue_ids_and_ixs(
            residue_ids, residue_ixs, "text"
        )
        assert residue_ids == residue_ids_formatted
        assert residue_ixs == residue_ixs_formatted
