"""
Useful utilities used across the library
"""

import os
import shutil
import tempfile
import contextlib


@contextlib.contextmanager
def enter_temp_directory(remove=True):
    """Create and enter a temporary directory; used as context manager."""
    temp_dir = tempfile.mkdtemp()
    cwd = os.getcwd()
    os.chdir(temp_dir)
    yield cwd, temp_dir
    os.chdir(cwd)
    if remove:
        shutil.rmtree(temp_dir)
