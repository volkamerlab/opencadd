"""
Tests for opencadd.structure.superposition.engines.mmligner
"""

import os
import pytest
from opencadd.structure.superposition.api import Structure
from opencadd.structure.superposition.engines.mmligner import MMLignerAligner


def test_mmligner_instantiation():
    aligner = MMLignerAligner()


@pytest.mark.skipif(
    "GITHUB_ACTIONS" in os.environ, reason="FIXME: MMLigner in conda-forge is not yet patched"
)
def test_mmligner_alignment():
    aligner = MMLignerAligner()
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]
    selections = []
    result = aligner.calculate(structures, selections)
    # Check API compliance
    assert "superposed" in result
    assert "scores" in result
    assert "rmsd" in result["scores"]
    assert "metadata" in result

    # Check RMSD value for these PDBs
    assert pytest.approx(result["scores"]["rmsd"], 2.279)


def test_mmligner_alignment_selections():
    aligner = MMLignerAligner()
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in ["4u3y", "4u40"]]
    user_select = ["backbone and name CA and segid A", "backbone and name CA and segid A"]
    selections = [
        structures[0].select_atoms(f"{user_select[0]}"),
        structures[1].select_atoms(f"{user_select[0]}"),
    ]
    result = aligner.calculate(structures, selections)
    # Check API compliance
    assert "superposed" in result
    assert "scores" in result
    assert "rmsd" in result["scores"]
    assert "metadata" in result

    # Check RMSD value for these PDBs
    assert pytest.approx(result["scores"]["rmsd"], 2.279)
