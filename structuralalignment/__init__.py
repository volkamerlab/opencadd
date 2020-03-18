"""
StructuralAlignment
A Python library for molecular structural alignment and superposition
"""

# Add imports here
from .structuralalignment import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
