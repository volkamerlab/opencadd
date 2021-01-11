"""
Tests for opencadd.structure.pocket.anchor
"""

import pandas as pd
import pytest

from opencadd.structure.pocket import AnchorResidue


class TestAnchorResidue:
    """
    Test AnchorResidue class methods.
    """

    @pytest.mark.parametrize(
        "center, residue_id, residue_id_alternative, residue_ix, color",
        [([1, 1, 1], 1, 1, 1, "blue")],
    )
    def test_attributes_and_properties(
        self, center, residue_id, residue_id_alternative, residue_ix, color
    ):
        anchor_residue = AnchorResidue(
            center, residue_id, residue_id_alternative, residue_ix, color, "subpocket", "pocket"
        )
        if center:
            for i, j in zip(anchor_residue.center, center):
                assert pytest.approx(i, abs=1.0e-6) == j
        assert anchor_residue.residue_id == residue_id
        assert anchor_residue.residue_id_alternative == residue_id_alternative
        assert anchor_residue.residue_ix == residue_ix
        assert anchor_residue.color == color
        assert isinstance(anchor_residue.anchor_residue, pd.Series)
        assert anchor_residue.anchor_residue.index.to_list() == [
            "anchor_residue.color",
            "anchor_residue.id",
            "anchor_residue.id_alternative",
            "anchor_residue.ix",
            "anchor_residue.center",
        ]
