"""
Unit and regression test for the devmos_stuff package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import devmos_stuff


def test_devmos_stuff_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "devmos_stuff" in sys.modules
