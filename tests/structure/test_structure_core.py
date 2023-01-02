"""
Tests for opencadd.structure.code
"""

import pytest
from opencadd.structure.core import Structure


def test_fetch_pdb():
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]

    assert len(structures) == 2
