"""
Tests for opencadd.structure.pocket.subpocket
"""

import pandas as pd
import pytest

from opencadd.structure.pocket import Subpocket, AnchorResidue

ANCHOR_RESIDUE1 = AnchorResidue([1, 1, 1], "1", ["2"], "1", "blue")
ANCHOR_RESIDUE2 = AnchorResidue([2, 2, 2], "2", ["3"], "2", "blue")
ANCHOR_RESIDUE3 = AnchorResidue(None, "2", ["3"], "2", "blue")


class TestSubpocket:
    """
    Test Subpocket class methods.
    """

    @pytest.mark.parametrize(
        "anchor_residues, name, color, center",
        [
            ([ANCHOR_RESIDUE1, ANCHOR_RESIDUE2], "hinge", "blue", [1.5, 1.5, 1.5]),
            ([ANCHOR_RESIDUE1, ANCHOR_RESIDUE3], "hinge", "blue", None),
            ([ANCHOR_RESIDUE3, ANCHOR_RESIDUE3], "hinge", "blue", None),
        ],
    )
    def test_attributes_and_properties(self, anchor_residues, name, color, center):
        subpocket = Subpocket(anchor_residues, name, color)

        # Test attributes
        assert subpocket.name == name
        assert subpocket.color == color
        for i, j in zip(subpocket._anchor_residues, anchor_residues):
            assert i == j

        # Test subpocket center calculation
        # (testing attribute center and method _centroid at the same time)
        if center:
            for i, j in zip(subpocket.center, center):
                assert i == j
        else:
            assert subpocket.center == center

        # Test properties
        assert isinstance(subpocket.data, pd.Series)
        assert subpocket.data.index.to_list() == [
            "subpocket.name",
            "subpocket.color",
            "subpocket.center",
        ]
        assert isinstance(subpocket.data_anchor_residues, pd.DataFrame)
        assert subpocket.data_anchor_residues.columns.to_list() == [
            "subpocket.name",
            "anchor_residue.color",
            "anchor_residue.id",
            "anchor_residue.id_alternative",
            "anchor_residue.ix",
            "anchor_residue.center",
        ]
