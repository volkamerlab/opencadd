"""
opencadd
A Python library for structural cheminformatics
"""

from ._version import __version__
import logging

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
