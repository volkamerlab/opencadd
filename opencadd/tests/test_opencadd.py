"""
Unit and regression test for the opencadd package.
"""

# Import package, test suite, and other packages as needed
import opencadd  # pylint: disable=unused-import
import sys


def test_opencadd_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "opencadd" in sys.modules
