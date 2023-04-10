"""
Tests for opencadd.structure.pocket.region
"""

import pandas as pd
import pytest

from opencadd.structure.pocket import Region


class TestRegion:
    """
    Test TestRegion class methods.
    """

    @pytest.mark.parametrize(
        "name, residue_ids, residue_ixs, color",
        [("example", [1, 2], [10, 20], "blue")],
    )
    def test_attributes_and_properties(self, name, residue_ids, residue_ixs, color):
        region = Region(name, residue_ids, residue_ixs, color)
        assert region.name == name
        assert region.residue_ids == residue_ids
        assert region.residue_ixs == residue_ixs
        assert region.color == color
        assert isinstance(region.region, pd.DataFrame)
        assert region.region.columns.to_list() == [
            "region.name",
            "region.color",
            "residue.id",
            "residue.ix",
        ]
