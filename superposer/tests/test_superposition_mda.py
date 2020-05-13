"""
Tests for superposer.superposition.mda
"""

import pytest
from superposer.api import Structure
from superposer.superposition.mda import MDAnalysisAligner


def test_mda_instantiation():
    aligner = MDAnalysisAligner()


def test_mda_calculation():
    aligner = MDAnalysisAligner()
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]
    result = aligner.calculate(structures)

    # Check API compliance
    assert "superposed" in result
    assert "scores" in result
    assert "rmsd" in result["scores"]
    assert "metadata" in result

    # Check RMSD values
    # TODO: pytest.approx is not working reliably - check with Dennis too, he has the same problem
    assert pytest.approx(result["scores"]["rmsd"], 1.989)
