"""
Tests for opencadd.structure.pocket.core
"""

import pytest

from opencadd.structure.pocket import BasePocket


class TestBasePocket:
    """
    Test BasePocket class methods.
    """

    @pytest.mark.parametrize(
        "residue_ids, residue_ixs, residue_ids_formatted, residue_ixs_formatted",
        [
            ([1, 2, 3], None, [1, 2, 3], [None, None, None]),
            (["1", "2", "_", "_"], ["1", "2", "3", "4"], [1, 2, None, None], [1, 2, 3, 4]),
            (["1", "2", None, None], ["1", "2", "3", "4"], [1, 2, None, None], [1, 2, 3, 4]),
        ],
    )
    def test_format_residue_ids_and_ixs(
        self, residue_ids, residue_ixs, residue_ids_formatted, residue_ixs_formatted
    ):
        """
        Test formatting of user-input residue PDB IDs and residue indices.
        """

        base_pocket = BasePocket()
        residue_ids2, residue_ixs2 = base_pocket._format_residue_ids_and_ixs(
            residue_ids, residue_ixs, ""
        )
        assert residue_ids2 == residue_ids_formatted
        assert residue_ixs2 == residue_ixs_formatted

    @pytest.mark.parametrize(
        "residue_ids, residue_ixs",
        [
            ([1, 2, 3], [None, 2, 3]),  # Non-int-castable index (None)
            ([1, 2, 3], ["a", 2, 3]),  # Non-int-castable index
            ([1, 1, 2], None),  # Duplicated PDB IDs
            ([1, 2, 3], [1, 1, 2]),  # Duplicated indices
        ],
    )
    def test_format_residue_ids_and_ixs_raises(self, residue_ids, residue_ixs):
        """
        Test error handling when formatting user-input residue PDB IDs and
        residue indices.
        """

        with pytest.raises((ValueError, TypeError)):
            base_pocket = BasePocket()
            base_pocket._format_residue_ids_and_ixs(residue_ids, residue_ixs, "")

    @pytest.mark.parametrize(
        "residue_ids, residue_ixs, n_residues",
        [
            ([101, None], [1, 2], 2),
            ([101, None], [1, 2], 2),
            ([101, None], [None, None], 2),
        ],
    )
    def test_residues(self, residue_ids, residue_ixs, n_residues):
        """
        Test property residues.
        """
        base_pocket = BasePocket()
        base_pocket._residue_ids = residue_ids
        base_pocket._residue_ixs = residue_ixs
        assert base_pocket.residues.columns.to_list() == ["residue.id", "residue.ix"]
        assert (
            base_pocket.residues.index.to_list()
            == base_pocket.residues.reset_index().index.to_list()
        )
        assert base_pocket.residues.dtypes.to_list() == ["Int32", "Int32"]
        assert len(base_pocket.residues) == n_residues

    @pytest.mark.parametrize(
        "residue_ids, residue_ixs, residue_id, residue_ix",
        [
            ([101, None], [1, 2], 101, 1),  # Residue ID+index exist
            ([101, None], [1, 2], 102, None),  # Residue ID does not exist
            ([101, None], [None, None], 101, None),  # Residue ID maps to None
        ],
    )
    def test_residue_id2ix(self, residue_ids, residue_ixs, residue_id, residue_ix):
        """
        Test residue PDB ID to index mapping.
        """
        base_pocket = BasePocket()
        base_pocket._residue_ids = residue_ids
        base_pocket._residue_ixs = residue_ixs
        assert base_pocket._residue_id2ix(residue_id) == residue_ix

    @pytest.mark.parametrize(
        "residue_ids, residue_ixs, residue_ix, residue_id",
        [
            ([101, None], [1, 2], 1, 101),  # Residue index+ID exist
            ([101, None], [1, 2], 2, None),  # Residue index maps to None
            ([101, 102], [1, 2], 10, None),  # Residue index does not exist
        ],
    )
    def test_residue_ix2id(self, residue_ids, residue_ixs, residue_ix, residue_id):
        """
        Test residue index to PDB ID mapping.
        """
        base_pocket = BasePocket()
        base_pocket._residue_ids = residue_ids
        base_pocket._residue_ixs = residue_ixs
        assert base_pocket._residue_ix2id(residue_ix) == residue_id
