"""
Tests for opencadd.structure.pocket.utils
"""

import pytest

from opencadd.structure.pocket.utils import _format_residue_ids_and_labels


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