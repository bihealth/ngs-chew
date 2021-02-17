import pathlib

import pytest


@pytest.fixture
def path_tests():
    return pathlib.Path(__file__).parent
