"""
Unit and regression test for the superposer package.
"""

# Import package, test suite, and other packages as needed
import superposer  # pylint: disable=unused-import
import sys


def test_superposer_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "superposer" in sys.modules
