"""
Tests for opencadd.structure.pocket.core
"""

from pathlib import Path

import pandas as pd
import pytest

from opencadd.structure.pocket.core import Pocket
from opencadd.structure.pocket.subpocket import AnchorResidue
from opencadd.structure.pocket.utils import _format_residue_ids_and_labels

PATH_TEST_DATA = Path(__name__).parent / "opencadd/tests/data/pocket"


@pytest.mark.parametrize(
    "residue_ids, residue_labels, residue_ids_formatted, residue_labels_formatted",
    [
        ([1, 2, 3], [1, 2, 3], ["1", "2", "3"], ["1", "2", "3"]),
        ([1, 2], None, ["1", "2"], [None, None]),
    ],
)
def test_format_residue_ids_and_labels(
    residue_ids, residue_labels, residue_ids_formatted, residue_labels_formatted
):

    residue_ids, residue_labels = _format_residue_ids_and_labels(residue_ids, residue_labels)
    assert residue_ids == residue_ids_formatted
    assert residue_labels == residue_labels_formatted


class TestsPocket:
    """
    Test Pocket class methods.
    """

    @pytest.mark.parametrize(
        "filepath, name, residue_ids, residue_labels, centroid",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                "example kinase",
                [127, 128, 129],
                [1, 2, 3],
                [-1.178433, 23.859733, 45.091933],
            )
        ],
    )
    def test_from_file_text(self, filepath, name, residue_ids, residue_labels, centroid):

        # Class method from_file calls from_text
        pocket = Pocket.from_file(filepath, residue_ids, name, residue_labels)

        assert pocket.name == name
        assert pocket.residues.equals(
            pd.DataFrame(
                {
                    "residue.id": [str(i) for i in residue_ids],
                    "residue.label": [str(i) for i in residue_labels],
                }
            )
        )
        assert pytest.approx(pocket.centroid, centroid, abs=1.0e-6)

    @pytest.mark.parametrize(
        "filepath, pocket_residue_ids, residue_ids, n_ca_atoms",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                ["12"],  # TODO remove str
                1,
            ),  # KLIFS ID 7139
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                ["13"],  # TODO remove str
                1,
            ),
        ],
    )
    def test_ca_atoms(self, filepath, pocket_residue_ids, residue_ids, n_ca_atoms):

        pocket = Pocket.from_file(filepath, pocket_residue_ids)
        ca_atoms = pocket._ca_atoms(*residue_ids)
        assert len(ca_atoms) == n_ca_atoms

    @pytest.mark.parametrize(
        "filepath, residue_ids, anchor_residue_id, anchor_residue_id_alternative, anchor_residue_color, anchor_residue_center",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [12, 16, 85, 87],
                "12",  # TODO remove str
                None,
                "yellow",
                [8.1173, 16.3077, 51.996],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                "86",
                ["85", "87"],  # TODO remove str
                "yellow",
                [3.8648999, 26.4354, 42.65935],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                "17",
                ["16"],
                "yellow",
                [9.3363, 11.0014, 42.1329],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                "18",
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
            assert pytest.approx(anchor_residue.center[0], anchor_residue_center[0])

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, anchor_residue_ix, anchor_residue_id, anchor_residue_id_alternative, anchor_residue_color, anchor_residue_center",
        [
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",  # KLIFS ID 7139
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                "1",  # TODO remove str
                "12",  # TODO remove str
                None,
                "yellow",
                [8.1173, 16.3077, 51.996],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                "45",
                None,
                ["85", "87"],  # TODO remove str
                "yellow",
                [3.8648999, 26.4354, 42.65935],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                "6",
                None,
                ["16"],
                "yellow",
                [9.3363, 11.0014, 42.1329],
            ),
            (
                PATH_TEST_DATA / "NEK2_5m55_altB_chainA.pdb",
                [12, 16, 85, 87],
                [1, 5, 44, 46],
                "7",
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

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_ixs)
        anchor_residue = pocket._anchor_residue_by_residue_ix(
            anchor_residue_ix, anchor_residue_color
        )
        assert isinstance(anchor_residue, AnchorResidue)
        assert anchor_residue.residue_ix == anchor_residue_ix
        assert anchor_residue.residue_id == anchor_residue_id
        assert anchor_residue.residue_id_alternative == anchor_residue_id_alternative
        assert anchor_residue.color == anchor_residue_color
        if anchor_residue_center:
            assert pytest.approx(anchor_residue.center[0], anchor_residue_center[0])

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, residue_id, residue_ix",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                ["127", "128", "129"],  # TODO remove str
                ["1", "2", "3"],  # TODO remove str
                "127",
                "1",
            ),
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                ["127", "128", "129"],
                ["1", "2", "3"],
                1,
                None,
            ),
        ],
    )
    def test_residue_id2ix(self, filepath, residue_ids, residue_ixs, residue_id, residue_ix):

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_ixs)
        assert pocket._residue_id2ix(residue_id) == residue_ix

    @pytest.mark.parametrize(
        "filepath, residue_ids, residue_ixs, residue_ix, residue_id",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                ["127", "128", "129"],  # TODO remove str
                ["1", "2", "3"],  # TODO remove str
                "1",
                "127",
            ),
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                ["127", "128", "129"],  # TODO remove str
                ["1", "2", "3"],  # TODO remove str
                "4",
                None,
            ),
        ],
    )
    def test_residue_ix2id(self, filepath, residue_ids, residue_ixs, residue_ix, residue_id):

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_ixs)
        assert pocket._residue_ix2id(residue_ix) == residue_id

    @pytest.mark.parametrize(
        "filepath, name, residue_ids, color, residue_labels",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                "hinge",
                [127, 128, 129],
                "magenta",
                [46, 47, 48],
            )
        ],
    )
    def test_add_regions(self, filepath, name, residue_ids, color, residue_labels):

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_labels)
        pocket.add_region(name, residue_ids, color, residue_labels)

        n_region_residues = len(residue_ids)
        region = pd.DataFrame(
            {
                "region.name": [name] * n_region_residues,
                "region.color": [color] * n_region_residues,
                "residue.id": [str(i) for i in residue_ids],
                "residue.label": [str(i) for i in residue_labels],
            }
        )
        assert pocket.regions.equals(region)

    @pytest.mark.parametrize(
        "filepath, name, residue_ids, color, residue_labels, center",
        [
            (
                PATH_TEST_DATA / "AAK1_4wsq_altA_chainA_protein.mol2",
                "hinge_region",
                [73, 128, 193],
                "magenta",
                [16, 47, 80],
                [1.957633, 21.923767, 41.69003333],
            )
        ],
    )
    def test_add_subpocket(self, filepath, name, residue_ids, color, residue_labels, center):

        pocket = Pocket.from_file(filepath, residue_ids, "", residue_labels)
        pocket.add_subpocket(name, residue_ids, color, residue_labels)
        assert pocket.subpockets.columns.to_list() == [
            "subpocket.name",
            "subpocket.color",
            "subpocket.center",
        ]

        subpocket = pocket._subpockets[0]
        assert subpocket.name == name
        assert subpocket.color == color
        assert pytest.approx(subpocket.center, center, abs=1.0e-6)
