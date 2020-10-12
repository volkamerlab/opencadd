"""
Tests for opencadd.databases.klifs.api
"""

from pathlib import Path

from bravado.client import SwaggerClient
import pandas as pd
import pytest

from opencadd.databases.klifs import local, remote
from opencadd.databases.klifs.api import setup_local, setup_remote

PATH_TEST_DATA = Path(__name__).parent / "opencadd/tests/data/klifs"


def test_api_remote():
    """
    Test Session attributes for remote session.
    """
    session = setup_remote()
    assert isinstance(session.client, SwaggerClient)
    assert session.database is None
    assert isinstance(session.kinases, remote.Kinases)
    assert isinstance(session.ligands, remote.Ligands)
    assert isinstance(session.structures, remote.Structures)
    assert isinstance(session.bioactivities, remote.Bioactivities)
    assert isinstance(session.interactions, remote.Interactions)
    assert isinstance(session.coordinates, remote.Coordinates)


@pytest.mark.parametrize(
    "path_to_klifs_download, path_to_klifs_metadata",
    [(PATH_TEST_DATA, None), (PATH_TEST_DATA, PATH_TEST_DATA / "klifs_metadata.csv")],
)
def test_api_local(path_to_klifs_download, path_to_klifs_metadata):
    """
    Test Session attributes for local session.

    Parameters
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder with KLIFS download files.
    path_to_klifs_metadata : pathlib.Path or str
        Path to KLIFS metadata file (default is None).
        Set this parameter, if you have initialized a local session before and therefore
        already have a KLIFS metadata file.
        You could pass here a filtered version of this KLIFS metadata file.
    """
    session = setup_local(path_to_klifs_download, path_to_klifs_metadata)

    assert session.client is None
    assert isinstance(session.database, pd.DataFrame)
    assert isinstance(session.kinases, local.Kinases)
    assert isinstance(session.ligands, local.Ligands)
    assert isinstance(session.structures, local.Structures)
    assert isinstance(session.bioactivities, local.Bioactivities)
    assert isinstance(session.interactions, local.Interactions)
    assert isinstance(session.coordinates, local.Coordinates)
