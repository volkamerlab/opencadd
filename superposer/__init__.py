"""
superposer
A Python library for molecular structural alignment and superposition
"""

# Handle versioneer
from .api import align, METHODS
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
