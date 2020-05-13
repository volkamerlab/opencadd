"""
Tests for superposer.api
"""

import pytest
from superposer.api import align, METHODS, Structure


@pytest.mark.parametrize("method", list(METHODS.values()))
def test_api(method):
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]
    results = align(structures, method=method)
    assert len(results) == len(structures) - 1
