"""
Tests for opencadd.structure.subpockets.core
"""

import pandas as pd
import pytest

from opencadd.structure.subpockets.core import Base, AnchorResidue

DATAFRAME = pd.DataFrame(
    {
        "residue.pdb_id": ["1", "2", "3", "7", "8", "9", "11"],
        "atom.name": ["CA", "CA", "CA", "CA", "CA", "CA", "CA"],
        "atom.x": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
        "atom.y": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
        "atom.z": [1.0, 2.0, 3.0, 4.0, 20.0, 30.0, 40.0],
    }
)


class TestsBase:
    """
    Test Base class methods.
    """

    @pytest.mark.parametrize("color_name, color_rgb", [("red", (1, 0, 0))])
    def test_format_color(self, color_name, color_rgb):
        base = Base()
        name, rgb = base._format_color(color_name)
        assert name == color_name
        assert rgb == color_rgb

    @pytest.mark.parametrize("color_name", ["xxx"])
    def test_format_color_raise(self, color_name):

        with pytest.raises(ValueError):
            base = Base()
            base._format_color(color_name)

    @pytest.mark.parametrize(
        "residue_pdb_ids, residue_labels, residue_pdb_ids_formatted, residue_labels_formatted",
        [
            ([1, 2, 3], [1, 2, 3], ["1", "2", "3"], ["1", "2", "3"]),
            ([1, 2], None, ["1", "2"], [None, None]),
        ],
    )
    def test_format_residue_pdb_ids_and_labels(
        self, residue_pdb_ids, residue_labels, residue_pdb_ids_formatted, residue_labels_formatted
    ):

        base = Base()
        residue_pdb_ids, residue_labels = base._format_residue_pdb_ids_and_labels(
            residue_pdb_ids, residue_labels
        )
        assert residue_pdb_ids == residue_pdb_ids_formatted
        assert residue_labels == residue_labels_formatted


class TestsAnchorResidue:
    """
    Test AnchorResidue class methods.
    """

    @pytest.mark.parametrize(
        "residue_pdb_id, residue_pdb_id_alternative, residue_center",
        [
            (1, None, [1, 1, 1]),
            (10, ["9", "11"], [35, 35, 35]),
            (4, ["3"], [3, 3, 3]),
            (6, ["7"], [4, 4, 4]),
            (5, None, None),
        ],
    )
    def test_from_dataframe(self, residue_pdb_id, residue_pdb_id_alternative, residue_center):
        """
        TODO
        """

        residue = AnchorResidue()
        residue.from_dataframe(DATAFRAME, residue_pdb_id)

        assert residue.pdb_id_alternative == residue_pdb_id_alternative
        if residue_center:
            assert pytest.approx(residue.center, residue_center, abs=1.0e-6)
