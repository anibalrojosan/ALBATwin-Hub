import pytest


@pytest.fixture
def sample_fixture():
    """Sample fixture for testing configuration."""
    return {"status": "ok"}
