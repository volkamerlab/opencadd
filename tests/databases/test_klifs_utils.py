"""
Tests for opencadd.databases.klifs.utils
"""

import pytest

from opencadd.db.klifs.utils import (
    metadata_to_filepath,
    filepath_to_metadata,
)


@pytest.mark.parametrize(
    "path_to_klifs_download, species, kinase_name, structure_pdb, structure_alternate_model, structure_chain, entity, extension, in_dir, filepath_template",
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
    extension,
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
        extension,
        in_dir,
    )
    assert str(filepath) == filepath_template


@pytest.mark.parametrize(
    "filepath, metadata_template",
    [
        (
            "/path/to/HUMAN/BRAF/6uuo/pocket.mol2",
            {
                "species": "Human",
                "kinase_name": "BRAF",
                "structure_pdb": "6uuo",
                "structure_alternate_model": None,
                "structure_chain": None,
                "entity": "pocket",
                "extension": "mol2",
            },
        ),
        (
            "/path/to/HUMAN/BRAF/6uuo_chainA/pocket.mol2",
            {
                "species": "Human",
                "kinase_name": "BRAF",
                "structure_pdb": "6uuo",
                "structure_alternate_model": None,
                "structure_chain": "A",
                "entity": "pocket",
                "extension": "mol2",
            },
        ),
        (
            "/path/to/HUMAN/BRAF/6uuo_altB_chainA/pocket.mol2",
            {
                "species": "Human",
                "kinase_name": "BRAF",
                "structure_pdb": "6uuo",
                "structure_alternate_model": "B",
                "structure_chain": "A",
                "entity": "pocket",
                "extension": "mol2",
            },
        ),
    ],
)
def test_filepath_to_metadata(filepath, metadata_template):
    """
    Test filepath to metadata conversion.
    """
    metadata = filepath_to_metadata(filepath)
    assert metadata == metadata_template
