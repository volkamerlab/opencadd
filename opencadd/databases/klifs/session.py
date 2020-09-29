"""
opencadd.databases.klifs.session

Defines the local/remote sessions.
"""

from pathlib import Path

import pandas as pd

from . import remote
from . import local
from .utils import KLIFS_CLIENT


class Session:
    """
    Class to set up remote or local session. Class attributes are None upon initialization but will
    be updated when remote or local session is called.

    Attributes
    ----------
    path_to_klifs_download : pathlib.Path
        Path to folder with KLIFS download files.
    client : None or bravado.client.SwaggerClient
        KLIFS client (set if session type is remote).
    database : None or pandas.DataFrame
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

    def __init__(self, backend, path_to_klifs_download=None, path_to_klifs_metadata=None):
        """
        TODO
        """

        # Not None in local only
        self.path_to_klifs_download = None
        self.database = None

        # Not None in remote only
        self.client = None

        # Not None in local and remote
        self.kinases = None
        self.ligands = None
        self.structures = None
        self.bioactivities = None
        self.interactions = None
        self.pockets = None
        self.coordinates = None

        # If path to a KLIFS download is set, initialize a local session
        if path_to_klifs_download:
            self.path_to_klifs_download = Path(path_to_klifs_download)

            if path_to_klifs_metadata:
                self.database = pd.read_csv(path_to_klifs_metadata)
            else:
                # TODO make simpler!
                session_initializer = local.SessionInitializer(self.path_to_klifs_download)
                self.database = session_initializer.klifs_metadata

        else:
            self.client = KLIFS_CLIENT

        # Set attributes using a backend (i.e. the local or remote module)
        self.kinases = backend.Kinases(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )
        self.ligands = backend.Ligands(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )
        self.structures = backend.Structures(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )
        self.bioactivities = backend.Bioactivities(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )
        self.interactions = backend.Interactions(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )
        self.pockets = backend.Pockets(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )
        self.coordinates = backend.Coordinates(
            client=self.client,
            database=self.database,
            path_to_klifs_download=self.path_to_klifs_download,
        )

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
        """

        return cls(local, path_to_klifs_download, path_to_klifs_metadata)

    @classmethod
    def from_remote(cls):
        """
        Set up remote session using the KLIFS swagger client.
        """

        return cls(remote)

