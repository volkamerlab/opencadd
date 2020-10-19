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


class TestsAnchorResidue:
    """
    Test AnchorResidue class methods.
    """

    @pytest.mark.parametrize(
        "residue_id, residue_id_alternative, residue_center",
        [
            (1, None, [1, 1, 1]),
            (10, ["9", "11"], [35, 35, 35]),
            (4, ["3"], [3, 3, 3]),
            (6, ["7"], [4, 4, 4]),
            (5, None, None),
        ],
    )
    def test_from_dataframe(self, residue_id, residue_id_alternative, residue_center):
        """
        Test anchor residue cases.
        """

        dataframe = pd.DataFrame(
            {
                "residue.id": ["1", "2", "3", "7", "8", "9", "11"],
                "atom.name": ["CA", "CA", "CA", "CA", "CA", "CA", "CA"],
                "atom.x": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
                "atom.y": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
                "atom.z": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
            }
        )

        residue = AnchorResidue.from_dataframe(dataframe, residue_id)

        assert residue.id_alternative == residue_id_alternative
        if residue_center:
            assert pytest.approx(residue.center, residue_center, abs=1.0e-6)


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
    def test_from_file(self, filepath, name, residue_ids, residue_labels, centroid):

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

        print(residue_ids)
        print(residue_labels)
        pocket = Pocket.from_file(filepath, residue_ids, "", residue_labels)

        pocket.add_subpocket(name, residue_ids, color, residue_labels)
        subpocket = pocket._subpockets[0]

        assert subpocket.name == name
        assert subpocket.color == color
        assert pytest.approx(subpocket.center, center, abs=1.0e-6)
        assert pocket.subpockets.columns.to_list() == [
            "subpocket.name",
            "subpocket.color",
            "subpocket.center",
        ]


# TODO Test KlifsPocket
