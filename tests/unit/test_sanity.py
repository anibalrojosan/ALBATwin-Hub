def test_sanity():
    """Verify that the testing environment is working correctly."""
    assert 1 + 1 == 2


def test_package_import():
    """Verify that the package can be imported."""
    import bioprocess_twin

    assert bioprocess_twin.__version__ == "0.1.0"
