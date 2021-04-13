"""
This module contains helper and utility functions.
"""
from pathlib import Path

__all__ = ["data_path"]


def data_path(fn):
    """Leads to files saved in the data folder

    Parameters
    ----------
    fn: str
        The whole filename.

    Returns
    -------
    The path of the file in the current working system.
    """
    return Path(__file__).parent / "data" / fn
