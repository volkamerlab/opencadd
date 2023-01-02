"""
Tests for opencadd.databases.klifs.exceptions
"""

from pathlib import Path

import pytest

from opencadd.api.klifs import setup_local
from opencadd.api.klifs.exceptions import (
    KlifsPocketIncompleteError,
    KlifsPocketUnequalSequenceStructure,
)

PATH_TEST_DATA = Path(__name__).parent / "opencadd/tests/data/klifs"

# Set local session
LOCAL = setup_local(PATH_TEST_DATA)


@pytest.mark.parametrize(
    "klifs_session, structure_klifs_id",
    [(LOCAL, 12347)],
)
def test_not_raise_error(klifs_session, structure_klifs_id):
    """
    Baseline test: Do not raise any error.
    """
    klifs_session.pockets.by_structure_klifs_id(structure_klifs_id, "pdb")


@pytest.mark.parametrize(
    "klifs_session, structure_klifs_id",
    [(LOCAL, 13623)],
)
def test_raise_klifs_pocket_incomplete_error(klifs_session, structure_klifs_id):
    """
    Example structure that has an incomplete pocket: Raise custom error!
    """
    with pytest.raises(KlifsPocketIncompleteError):
        klifs_session.pockets.by_structure_klifs_id(structure_klifs_id, "pdb")


@pytest.mark.parametrize(
    "klifs_session, structure_klifs_id",
    [(LOCAL, 1243)],
)
def test_raise_klifs_pocket_unequal_sequence_structure(klifs_session, structure_klifs_id):
    """
    Example structure has different pocket lengths in sequence and structure: Raise custom error!
    """
    with pytest.raises(KlifsPocketUnequalSequenceStructure):
        klifs_session.pockets.by_structure_klifs_id(structure_klifs_id, "pdb")
