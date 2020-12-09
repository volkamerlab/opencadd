"""
Tests for opencadd.structure.pocket.utils
"""

import pytest

from opencadd.structure.pocket.utils import _format_residue_ids_and_ixs


@pytest.mark.parametrize(
    "residue_ids, residue_ixs, residue_ids_formatted, residue_ixs_formatted",
    [
        ([1, 2, 3], [1, 2, 3], ["1", "2", "3"], ["1", "2", "3"]),
        ([1, 2], None, ["1", "2"], [None, None]),
    ],
)
def test_format_residue_ids_and_ixs(
    residue_ids, residue_ixs, residue_ids_formatted, residue_ixs_formatted
):

    residue_ids, residue_ixs = _format_residue_ids_and_ixs(residue_ids, residue_ixs)
    assert residue_ids == residue_ids_formatted
    assert residue_ixs == residue_ixs_formatted
