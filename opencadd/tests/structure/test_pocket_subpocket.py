"""
Tests for opencadd.structure.pocket.subpocket
"""

import pandas as pd
import pytest

from opencadd.structure.pocket import Subpocket, AnchorResidue

ANCHOR_RESIDUE1 = AnchorResidue([1, 1, 1], "1", ["2"], "1", "blue")
ANCHOR_RESIDUE2 = AnchorResidue([2, 2, 2], "2", ["3"], "2", "blue")


class TestSubpocket:
    """
    Test Subpocket class methods.
    """

    @pytest.mark.parametrize(
        "anchor_residues, name, color, center",
        [([ANCHOR_RESIDUE1, ANCHOR_RESIDUE2], "hinge", "blue", [1.5, 1.5, 1.5])],
    )
    def test_from_anchor_residue(self, anchor_residues, name, color, center):
        subpocket = Subpocket.from_anchor_residues(anchor_residues, name, color)

        # Test attributes
        assert subpocket.name == name
        assert subpocket.color == color
        for i, j in zip(subpocket.center, center):
            assert i == j
        for i, j in zip(subpocket._anchor_residues, anchor_residues):
            assert i == j

        # Test properties
        assert isinstance(subpocket.data, pd.Series)
        assert subpocket.data.index.to_list() == [
            "subpocket.name",
            "subpocket.color",
            "subpocket.center",
        ]
        assert isinstance(subpocket.anchor_residues, pd.DataFrame)
        assert subpocket.anchor_residues.columns.to_list() == [
            "subpocket.name",
            "subpocket.color",
            "anchor_residue.id",
            "anchor_residue.id_alternative",
            "anchor_residue.ix",
            "anchor_residue.center",
        ]

        # Test method
        for i, j in zip(subpocket._centroid(), center):
            assert i == j
