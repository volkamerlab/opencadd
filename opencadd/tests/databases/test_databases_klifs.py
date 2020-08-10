"""
Unit and regression test for the opencadd.databases.klifs module.
"""

# Import package, test suite, and other packages as needed
import opencadd.databases.klifs
import pytest
import sys


def test_databases_klifs_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "opencadd.databases.klifs" in sys.modules
