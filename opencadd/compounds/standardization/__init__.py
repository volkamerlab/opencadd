"""
standardization
A tool to standardize compounds.
"""

# Add imports here
from .standardization import *
from .assign_stereochemistry import *
from .convert_format import *
from .detect_inorganic import *
from .disconnect_metals import *
from .handle_charges import *
from .handle_fragments import *
from .handle_hydrogens import *
from .normalize_molecules import *
from .remove_salts import *
from .utils import *
from .sanitize_molecules import *
from .validate_molecules import *


# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
