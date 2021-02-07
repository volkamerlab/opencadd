"""
standardization
A tool to standardize compounds.
"""

# Add imports here
from opencadd.compounds.standardization.standardization import *
from opencadd.compounds.standardization.assign_stereochemistry import *
from opencadd.compounds.standardization.convert_format import *
from opencadd.compounds.standardization.detect_inorganic import *
from opencadd.compounds.standardization.disconnect_metals import *
from opencadd.compounds.standardization.handle_charges import *
from opencadd.compounds.standardization.handle_fragments import *
from opencadd.compounds.standardization.handle_hydrogens import *
from opencadd.compounds.standardization.normalize_molecules import *
from opencadd.compounds.standardization.remove_salts import *
from opencadd.compounds.standardization.utils import *
from opencadd.compounds.standardization.sanitize_molecules import *
from opencadd.compounds.standardization.validate_molecules import *


# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
