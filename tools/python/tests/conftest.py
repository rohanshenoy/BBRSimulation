"""Shared fixtures for the bbrsim tests.

The directory has no __init__.py, so the tests import whichever bbrsim is on
PYTHONPATH (under the runner: the installed copy). Run them from / or the
repository root, never from tools/python, where the source package would
shadow it.
"""
import os
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def repo_root():
    """The source tree holding tools/python/tests (validation/ and data/ beside it)."""
    return Path(__file__).resolve().parents[3]


@pytest.fixture(scope="session")
def data_root():
    """$BBRSIMDATA (the runtime data root); skips the test when it is not set."""
    d = os.environ.get("BBRSIMDATA")
    if not d:
        pytest.skip("BBRSIMDATA is not set")
    return d
