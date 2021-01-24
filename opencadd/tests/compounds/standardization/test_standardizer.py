"""
Unit and regression test for the standardizer package.
"""

# Import package, test suite, and other packages as needed

import sys
import pytest
import opencadd.compounds.standardization


def test_standardizer_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "standardizer" in sys.modules
