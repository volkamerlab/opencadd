"""
api.py

Defines the opencadd.databases.klifs API (local and remote).
"""
import logging

import pandas as pd

from . import remote
from . import local
from .utils import KLIFS_CLIENT


# Set logger
_logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

# Needed to show all columns in Jupyter notebooks
pd.set_option("display.max_columns", 100)


def setup_remote():
    """
    Set up remote session to work with KLIFS data remotely.

    Returns
    -------
    opencadd.databases.klifs.api.Session
        Remote session.
    """
    _logger.info(f"Set up remote session...")
    session = Session()
    session.from_remote()
    _logger.info(f"Remote session is ready!")
    return session


def setup_local(path_to_klifs_download):
    """
    Set up local session to work with KLIFS data locally.

    Parameter
    ---------
    path_to_klifs_download : pathlib.Path or str
        Path to folder with KLIFS download files.
    
    Returns
    -------
    opencadd.databases.klifs.api.Session
        Remote session.
    """
    _logger.info(f"Set up local session...")
    session = Session()
    session.from_local(path_to_klifs_download)
    _logger.info(f"Local session is ready!")
    return session


class Session:
    """
    Class to set up remote or local session. Class attributes are None upon initialization but will
    be updated when remote or local session is called.

    Attributes
    ----------
    session_type : None or str
        Session type, i.e. remote or local. 
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

    def __init__(self):

        self.session_type = None
        self.client = None
        self.database = None
        self.kinases = None
        self.ligands = None
        self.structures = None
        self.bioactivities = None
        self.interactions = None
        self.pockets = None
        self.coordinates = None

    def from_local(self, path_to_klifs_download):
        """
        Set up local session by initializing a local metadata database from KLIFS download.

        Parameters
        ----------
        path_to_klifs_download : pathlib.Path or str
            Path to folder with KLIFS download files.
        """

        # Session type
        self.session_type = "local"

        # Get database
        session_initializer = local.SessionInitializer(path_to_klifs_download)
        self.database = session_initializer.klifs_metadata

        # Initialize classes
        self.kinases = local.Kinases(self.database)
        self.ligands = local.Ligands(self.database)
        self.structures = local.Structures(self.database)
        self.bioactivities = local.Bioactivities(self.client)
        self.interactions = local.Interactions(self.database)
        self.pockets = local.Pockets(self.database)
        self.coordinates = local.Coordinates(self.database)

    def from_remote(self):
        """
        Set up remote session using the KLIFS swagger client.
        """

        # Session type
        self.session_type = "remote"

        # Get client
        self.client = KLIFS_CLIENT

        # Initialize classes
        self.kinases = remote.Kinases(self.client)
        self.ligands = remote.Ligands(self.client)
        self.structures = remote.Structures(self.client)
        self.bioactivities = remote.Bioactivities(self.client)
        self.interactions = remote.Interactions(self.client)
        self.pockets = remote.Pockets(self.client)
        self.coordinates = remote.Coordinates(self.client)
