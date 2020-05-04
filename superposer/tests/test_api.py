"""
Tests for superposer.api
"""

import pytest
import atomium
from superposer.api import align, METHODS


@pytest.mark.parametrize("method", list(METHODS.values()))
def test_api(method):
    structures = [atomium.fetch(pdb_id).model for pdb_id in ["4u3y", "4u40"]]
    results = align(structures, method=method)
    assert len(results) == len(structures) - 1
