"""
Tests for opencadd.databases.klifs.utils
"""

import pytest

from opencadd.databases.klifs_new.utils import (
    metadata_to_filepath,
    filepath_to_metadata,
)


@pytest.mark.parametrize(
    "path_to_klifs_download, species, kinase_name, structure_pdb, structure_alternate_model, structure_chain, entity, format, in_dir, filepath_template",
    [
        (
            "/path/to",
            "human",
            "BRAF",
            "6uuo",
            None,
            None,
            "complex",
            "mol2",
            True,
            "/path/to/HUMAN/BRAF/6uuo/complex.mol2",
        ),
        (
            "/path/to",
            "human",
            "BRAF",
            "6uuo",
            None,
            "A",
            "complex",
            "mol2",
            True,
            "/path/to/HUMAN/BRAF/6uuo_chainA/complex.mol2",
        ),
        (
            "/path/to",
            "human",
            "BRAF",
            "6uuo",
            "B",
            "A",
            "complex",
            "mol2",
            True,
            "/path/to/HUMAN/BRAF/6uuo_altB_chainA/complex.mol2",
        ),
        (
            "/path/to",
            "human",
            "BRAF",
            "6uuo",
            "B",
            "A",
            "complex",
            "mol2",
            False,
            "/path/to/HUMAN_BRAF_6uuo_altB_chainA_complex.mol2",
        ),
    ],
)
def test_metadata_to_filepath(
    path_to_klifs_download,
    species,
    kinase_name,
    structure_pdb,
    structure_alternate_model,
    structure_chain,
    entity,
    format,
    in_dir,
    filepath_template,
):
    """
    Test metadata to filepath conversion.
    """

    filepath = metadata_to_filepath(
        path_to_klifs_download,
        species,
        kinase_name,
        structure_pdb,
        structure_alternate_model,
        structure_chain,
        entity,
        format,
        in_dir,
    )
    assert str(filepath) == filepath_template

