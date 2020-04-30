"""
Tests for structuralalignment.superposition.mmligner
"""

import pytest
import atomium
from structuralalignment.superposition.mmligner import MMLignerAligner


def test_mmligner_instantiation():
    aligner = MMLignerAligner()


def test_mmligner_alignment():
    aligner = MMLignerAligner()
    structures = [atomium.fetch(pdb_id).model for pdb_id in ["4u3y", "4u40"]]
    result = aligner.calculate(structures)
    # Check API compliance
    assert "superposed" in result
    assert "scores" in result
    assert "rmsd" in result["scores"]
    assert "metadata" in result

    # Check RMSD value for these PDBs
    assert pytest.approx(result["scores"]["rmsd"], 2.279)


# @pytest.mark.parametrize("structure1, structure2, structure3", ["ARN", "DBC", "EQZGHILKM"])
# def test_calculate_invalid_inputs(structure1, structure2, structure3):
#     mmligner = MMLignerAligner()

#     with pytest.raises(TypeError):
#         mmligner.calculate([[structure1], [None]])

#     with pytest.raises(TypeError):
#         mmligner.calculate([None, structure2])

#     with pytest.raises(TypeError):
#         mmligner.calculate([structure1])

#     with pytest.raises(TypeError):
#         mmligner.calculate([structure1, structure2, structure3])
