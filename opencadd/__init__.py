"""
opencadd
A Python library for structural cheminformatics
"""

import logging

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

_logger = logging.getLogger(__name__)

_logger.info(
    """
    Please note the following development stages of opencadd's submodules:

    Pre-release
    -----------
    - opencadd.databases.klifs
    - opencadd.structure.pocket

    Work-in-progress
    ----------------
    - opencadd.io
    - opencadd.structure.superposition
    """
)


from .protein import Protein