"""
opencadd.databases.klifs.api

Defines the opencadd.databases.klifs API (local and remote).
"""

import logging

import pandas as pd

from .session import Session

# Set logger
_logger = logging.getLogger(__name__)

# Add warning on truncated DataFrames in Jupyter notebooks
if pd.get_option("display.max_columns") < 50:
    _logger.info(
        "If you want to see an non-truncated version of the DataFrames in this module, "
        "use `pd.set_option('display.max_columns', 50)` in your notebook."
    )


def setup_remote():
    """
    Set up remote session to work with KLIFS data remotely.

    Returns
    -------
    opencadd.databases.klifs.session.Session
        Remote session.
    """
    _logger.info("Set up remote session...")
    session = Session.from_remote()
    _logger.info("Remote session is ready!")
    return session


def setup_local(path_to_klifs_download, path_to_klifs_metadata=None):
    """
    Set up local session to work with KLIFS data locally.

    Parameter
    ---------
    path_to_klifs_download : pathlib.Path or str
        Path to folder with KLIFS download files.
    path_to_klifs_metadata : pathlib.Path or str
        Path to KLIFS metadata file (default is None).
        Set this parameter, if you have initialized a local session before and therefore
        already have a KLIFS metadata file.
        You could pass here a filtered version of this KLIFS metadata file.

    Returns
    -------
    opencadd.databases.klifs.session.Session
        Remote session.
    """
    _logger.info("Set up local session...")
    session = Session.from_local(path_to_klifs_download, path_to_klifs_metadata)
    _logger.info("Local session is ready!")
    return session
