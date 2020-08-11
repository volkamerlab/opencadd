"""
api.py

Defines API.
"""

from bravado.client import SwaggerClient

from . import remote
from . import local


KLIFS_API_DEFINITIONS = "http://klifs.vu-compmedchem.nl/swagger/swagger.json"
KLIFS_CLIENT = SwaggerClient.from_url(
    KLIFS_API_DEFINITIONS, config={"validate_responses": False}
)


class Session:
    """Class for remote or local session"""

    def __init__(self, path_to_klifs_download=None):

        if path_to_klifs_download:
            # Session type
            self.session_type = "local"

            # Get database
            session_initializer = local.SessionInitializer(path_to_klifs_download)
            self.database = session_initializer.klifs_metadata

            # Initialize classes
            self.kinases = local.Kinases(self.database)
            self.ligands = local.Ligands(self.database)
            self.structures = local.Structures(self.database)
            self.interactions = local.Interactions(self.database)
            self.coordinates = local.Coordinates(self.database)

        else:
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

