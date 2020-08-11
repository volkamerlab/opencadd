"""
api.py

Defines opencadd.databases.klifs API (local and remote).
"""

from bravado.client import SwaggerClient

from . import remote
from . import local


KLIFS_API_DEFINITIONS = "http://klifs.vu-compmedchem.nl/swagger/swagger.json"
KLIFS_CLIENT = SwaggerClient.from_url(
    KLIFS_API_DEFINITIONS, config={"validate_responses": False}
)


class Session:
    """
    Class to set up remote or local session.
    Class attributes are None upon initialization but will be updated when remote or local session is called, 
    using from_remote() or from_local(path_to_klifs_download).

    Attributes
    ----------
    session_type : None or str
        Session type, i.e. remote or local. 
    client : None or bravado.client.SwaggerClient
        KLIFS client (set if session type is remote).
    database : None or pandas.DataFrame
        KLIFS metadata (set if session type is local).
    kinases : None or remote.Kinases or local.Kinases
        Kinases object for kinases requests.
    ligands : None or remote.Ligands or local.Ligands
        Ligands object for ligands requests.
    structures: None or remote.Structures or local.Structures
        Structures object for structures requests.
    bioactivities : None or remote.Bioactivities or local.Bioactivities
        Bioactivities object for bioactivity requests.
    interactions : None or remote.Interactions or local.Interactions
        Interactions object for interactions requests.
    coordinates : None or remote.Coordinates or local.Coordinates
        Coordinates object for coordinates requests.
    """

    def __init__(self):

        self.session_type = None
        self.client = None
        self.database = None
        self.kinases = None
        self.kinases = None
        self.ligands = None
        self.structures = None
        self.bioactivities = None
        self.interactions = None
        self.coordinates = None

    def from_local(self, path_to_klifs_download):
        """
        Set up local session by initializing a local metadata database from KLIFS download.
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
        self.coordinates = remote.Coordinates(self.client)
