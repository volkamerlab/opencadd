"""
Tests for opencadd.databases.klifs.api
"""

from pathlib import Path

from bravado.client import SwaggerClient
import pandas as pd
import pytest

from opencadd.databases.klifs_new import local
from opencadd.databases.klifs_new import remote
from opencadd.databases.klifs_new.api import setup_local, setup_remote

PATH_TEST_DATA = Path(__file__).parent / "data"


def test_api_remote():
    """
    Test Session attributes for remote session.
    """
    session = setup_remote()
    assert session.session_type == "remote"
    assert isinstance(session.client, SwaggerClient)
    assert session.database is None
    assert isinstance(session.kinases, remote.Kinases)
    assert isinstance(session.ligands, remote.Ligands)
    assert isinstance(session.structures, remote.Structures)
    assert isinstance(session.bioactivities, remote.Bioactivities)
    assert isinstance(session.interactions, remote.Interactions)
    assert isinstance(session.coordinates, remote.Coordinates)


@pytest.mark.parametrize("path_to_klifs_download", [PATH_TEST_DATA / "KLIFS_download"])
def test_api_local(path_to_klifs_download):
    """
    Test Session attributes for local session.

    Parameters
    ----------
    path_to_klifs_download : pathlib.Path or str
        Path to folder with KLIFS download files.
    """
    session = setup_local(path_to_klifs_download)
    assert session.session_type == "local"
    assert session.client is None
    assert isinstance(session.database, pd.DataFrame)
    assert isinstance(session.kinases, local.Kinases)
    assert isinstance(session.ligands, local.Ligands)
    assert isinstance(session.structures, local.Structures)
    assert isinstance(session.bioactivities, local.Bioactivities)
    assert isinstance(session.interactions, local.Interactions)
    assert isinstance(session.coordinates, local.Coordinates)
