"""
Tests for opencadd.structure.pocket.core
"""

from pathlib import Path

import pandas as pd
import pytest

from opencadd.structure.pocket.core import Pocket, AnchorResidue

PATH_TEST_DATA = Path(__name__).parent / "opencadd/tests/data/pocket"


class TestsPocket:
    """
    Test Pocket class methods.
    """

    @pytest.mark.parametrize(
        "filepath, residue_ids, name, residue_ixs",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [44, 45],
                "example kinase",
                None,
            )
        ],
    )
    def test_from_file(self, filepath, residue_ids, name, residue_ixs):
        """
        Initialize class from file and test class attributes and properties.
        """

        pocket = Pocket.from_file(filepath, residue_ids, name, residue_ixs)

        # Test attributes
        assert pocket.name == name
        assert isinstance(pocket.data, pd.DataFrame)
        assert list(pocket.data["residue.id"].unique()) == [str(i) for i in residue_ids]
        assert isinstance(pocket._text, str)
        assert pocket._extension == filepath.suffix[1:]
        assert pocket._residue_ids == [str(i) for i in residue_ids]
        if residue_ixs:
            assert pocket._residue_ixs == [str(i) for i in residue_ixs]
        else:
            assert pocket._residue_ixs == [None] * len(residue_ids)

        # Test properties
        assert pocket.residues.columns.to_list() == ["residue.id", "residue.ix"]

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, region_name, region_residue_ids, region_color",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                "hinge",
                [127, 128, 129],
                "magenta",
            )
        ],
    )
    def test_add_regions(
        self, filepath, pocket_residue_ids, region_name, region_residue_ids, region_color
    ):
        """
        Test adding regions and associated regions property and attribute.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids)
        assert pocket.regions == None

        # Add region
        pocket.add_region(region_name, region_residue_ids, region_color)

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
        "filepath, pocket_residue_ids, subpocket_name, subpocket_residue_ids, subpocket_color, subpocket_center",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [73, 128, 193],
                "hinge_region",
                [73, 128, 193],
                "magenta",
                [1.957633, 21.923767, 41.69003333],
            )
        ],
    )
    def test_add_subpocket(
        self,
        filepath,
        pocket_residue_ids,
        subpocket_name,
        subpocket_residue_ids,
        subpocket_color,
        subpocket_center,
    ):
        """
        Test adding subpocket and associated subpockets property and attribute.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids)
        assert pocket.subpockets == None
        assert pocket.anchor_residues == None

        # Add subpocket
        pocket.add_subpocket(subpocket_name, subpocket_residue_ids, subpocket_color)

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

        # Test an example subpocket in the _subpockets attribute
        subpocket = pocket._subpockets[0]
        assert subpocket.name == subpocket_name
        assert subpocket.color == subpocket_color
        assert pytest.approx(subpocket.center, abs=1.0e-6) == subpocket_center

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

        pocket = Pocket.from_file(filepath, pocket_residue_ids, "", pocket_residue_ixs)
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
                ["85", "87"],
                "yellow",
                [3.8648999, 26.4354, 42.65935],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                17,
                ["16"],
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
        assert anchor_residue.residue_id == str(anchor_residue_id)
        assert anchor_residue.residue_id_alternative == anchor_residue_id_alternative
        assert anchor_residue.residue_ix == None
        assert anchor_residue.color == anchor_residue_color
        if anchor_residue_center:
            assert pytest.approx(anchor_residue.center[0], anchor_residue_center[0])

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, anchor_residue_ix, anchor_residue_id, anchor_residue_id_alternative, anchor_residue_color, anchor_residue_center",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                1,
                "12",
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
                ["85", "87"],
                "yellow",
                [3.8648999, 26.4354, 42.65935],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                6,
                None,
                ["16"],
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

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_ixs)
        anchor_residue = pocket._anchor_residue_by_residue_ix(
            anchor_residue_ix, anchor_residue_color
        )
        assert isinstance(anchor_residue, AnchorResidue)
        assert anchor_residue.residue_ix == str(anchor_residue_ix)
        if anchor_residue_id:
            assert anchor_residue.residue_id == str(anchor_residue_id)
        assert anchor_residue.residue_id_alternative == anchor_residue_id_alternative
        assert anchor_residue.color == anchor_residue_color
        if anchor_residue_center:
            assert pytest.approx(
                anchor_residue.center[0], anchor_residue_center[0]
            )  # TODO check this!

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, residue_id, residue_ix",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                [1, 2, 3],
                127,
                "1",
            ),
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                [1, 2, 3],
                1,
                None,
            ),
        ],
    )
    def test_residue_id2ix(self, filepath, residue_ids, residue_ixs, residue_id, residue_ix):
        """
        Test residue PDB ID to index mapping.
        """

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_ixs)
        assert pocket._residue_id2ix(residue_id) == residue_ix

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, residue_ix, residue_id",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                [1, 2, 3],
                1,
                "127",
            ),
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                [1, 2, 3],
                4,
                None,
            ),
        ],
    )
    def test_residue_ix2id(self, filepath, residue_ids, residue_ixs, residue_ix, residue_id):
        """
        Test residue index to PDB ID  mapping.
        """

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_ixs)
        assert pocket._residue_ix2id(residue_ix) == residue_id

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
        ca_atoms = pocket._ca_atoms(*residue_ids)
        assert len(ca_atoms) == n_ca_atoms
        _, ca_atoms_center = pocket._ca_atoms_center(*residue_ids)

        if center:
            for i, j in zip(ca_atoms_center, center):
                assert pytest.approx(i, abs=1.0e-3) == j

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, pocket_center",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                [127, 128, 129],
                [-1.178, 23.860, 45.092],
            ),
        ],
    )
    def test_center(self, filepath, pocket_residue_ids, pocket_center):
        """
        Test the pocket center calculation based on the pocket's residue PDB IDs.
        """

        pocket = Pocket.from_file(filepath, pocket_residue_ids)

        for i, j in zip(pocket.center, pocket_center):
            assert pytest.approx(i, abs=1.0e-3) == j
