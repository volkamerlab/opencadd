"""
Tests for opencadd.api
"""

import pytest
from opencadd.structure.superposition.api import align, METHODS, Structure


@pytest.mark.parametrize("method", list(METHODS.values()))
def test_api(method):
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]
    user_select = []
    results = align(structures, user_select, method=method)
    assert len(results) == len(structures) - 1


@pytest.mark.parametrize("method", list(METHODS.values()))
def test_api_selections(method):
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]
    user_select = ["backbone and name CA and segid A", "backbone and name CA and segid A"]
    results = align(structures, user_select, method=method)
    assert len(results) == len(structures) - 1
    assert len(user_select) == len(structures)
