"""
Unit and regression test for the standardization package.
"""

# Import package, test suite, and other packages as needed

import sys
import pytest
import opencadd.compounds.standardization


def test_standardization_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "opencadd" in sys.modules
