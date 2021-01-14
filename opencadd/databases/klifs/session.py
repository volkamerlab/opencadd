"""
opencadd.databases.klifs.session

Defines the local/remote sessions.
"""

from pathlib import Path

import pandas as pd

from . import remote
from . import local
from .remote import KLIFS_CLIENT


class Session:
    """
    Class to set up remote or local session.

    Attributes
    ----------
    _path_to_klifs_download : None or pathlib.Path
        Path to folder with KLIFS download files.
    _client : None or bravado.client.SwaggerClient
        KLIFS client (set if session type is remote).
    _database : None or pandas.DataFrame
        KLIFS metadata (set if session type is local).
    kinases : None or opencadd.databases.klifs.remote.Kinases/local.Kinases
        Kinases object for kinases requests.
    ligands : None or opencadd.databases.klifs.remote.Ligands/local.Ligands
        Ligands object for ligand requests.
    structures: None or opencadd.databases.klifs.remote.Structures/local.Structures
        Structures object for structure requests.
    bioactivities : None or opencadd.databases.klifs.remote.Bioactivities/local.Bioactivities
        Bioactivities object for bioactivity requests.
    interactions : None or opencadd.databases.klifs.remote.Interactions/local.Interactions
        Interactions object for interaction requests.
    coordinates : None or opencadd.databases.klifs.remote.Coordinates/local.Coordinates
        Coordinates object for coordinates requests.
    """

    def __init__(self):

        # Not None in local only
        self._path_to_klifs_download = None
        self._database = None

        # Not None in remote only
        self._client = None

        # Not None in local and remote
        self.kinases = None
        self.ligands = None
        self.structures = None
        self.bioactivities = None
        self.interactions = None
        self.pockets = None
        self.coordinates = None

    @classmethod
    def from_local(cls, path_to_klifs_download, path_to_klifs_metadata=None):
        """
        Set up local session by initializing or loading a local metadata database
        from a KLIFS download.

        Parameters
        ----------
        path_to_klifs_download : pathlib.Path or str
            Path to folder with KLIFS download files.
        path_to_klifs_metadata : pathlib.Path or str
            Path to KLIFS metadata file (default is None).
            Set this parameter, if you have initialized a local session before and therefore
            already have a KLIFS metadata file.
            You could pass here a filtered version of this KLIFS metadata file.

        Returns
        -------
        opencadd.databases.klifs.Session
            Local KLIFS session.
        """

        path_to_klifs_download = Path(path_to_klifs_download)
        if path_to_klifs_metadata:
            database = pd.read_csv(path_to_klifs_metadata)
        else:
            database = local._LocalDatabaseGenerator.from_files(path_to_klifs_download)

        session = cls()
        session._set_attributes(
            local, database=database, path_to_klifs_download=path_to_klifs_download
        )

        return session

    @classmethod
    def from_remote(cls):
        """
        Set up remote session using the KLIFS swagger client.

        Returns
        -------
        opencadd.databases.klifs.Session
            Remote KLIFS session.
        """

        client = KLIFS_CLIENT

        session = cls()
        session._set_attributes(remote, client=client)

        return session

    def _set_attributes(self, backend, path_to_klifs_download=None, database=None, client=None):
        """
        Set attributes using a backend (i.e. the local or remote module).

        Parameters
        ----------
        backend : opencadd.databases.klifs.local or opencadd.databases.klifs.remote
            Local or remote module.
        path_to_klifs_download : None or pathlib.Path
            Path to folder with KLIFS download files.
        client : None or bravado.client.SwaggerClient
            KLIFS client (set if session type is remote).
        database : None or pandas.DataFrame
            KLIFS metadata (set if session type is local).
        """

        self._client = client
        self._database = database
        self._path_to_klifs_download = path_to_klifs_download

        self.kinases = backend.Kinases(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
        self.ligands = backend.Ligands(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
        self.structures = backend.Structures(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
        self.bioactivities = backend.Bioactivities(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
        self.interactions = backend.Interactions(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
        self.pockets = backend.Pockets(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
        self.coordinates = backend.Coordinates(
            client=client,
            database=database,
            path_to_klifs_download=path_to_klifs_download,
        )
