"""
opencadd.databases.klifs.klifs_client
Utility functions to work with KLIFS data (remote)

KLIFS client.
"""

from bravado.client import SwaggerClient

KLIFS_API_DEFINITIONS = "http://klifs.vu-compmedchem.nl/swagger/swagger.json"
KLIFS_CLIENT = SwaggerClient.from_url(KLIFS_API_DEFINITIONS, config={"validate_responses": False})
