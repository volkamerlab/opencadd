"""
Test structuralalignment.superposition.theseus
"""

import pytest
import atomium


# TODO: Fill in the unit tests


def test_theseus_instantiation():
    from structuralalignment.superposition.theseus import TheseusAligner

    aligner = TheseusAligner()


def test_theseus_identical():
    pass


def test_theseus_different():
    from structuralalignment.superposition.theseus import TheseusAligner

    different_models = [atomium.fetch(pdb_id).model for pdb_id in ["6HG4", "6HG9"]]
    aligner = TheseusAligner()
    results = aligner.calculate(different_models, identical=False)
    # Check API compliance
    assert "superposed" in results
    assert "scores" in results
    assert "metadata" in results
    assert pytest.approx(results["scores"]["rmsd"], 1.54)


def test_parse_muscle():
    # TODO: Andrea Volkamer suggested we report the coverage
    # and this is part of the MUSCLE output
    pass


def test_parse_rmsd():
    pass


def test_parse_others():
    pass
