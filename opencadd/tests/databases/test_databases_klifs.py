"""
Unit and regression test for the opencadd.databases.klifs module.
"""

# Import package, test suite, and other packages as needed
import opencadd.databases.klifs
import pytest
import sys

def test_klifs_utils_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "klifs_utils" in sys.modules
