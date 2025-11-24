"""Unit tests for __version__.py."""

import diffpy.cmisas  # noqa


def test_package_version():
    """Ensure the package version is defined and not set to the initial
    placeholder."""
    assert hasattr(diffpy.cmisas, "__version__")
    assert diffpy.cmisas.__version__ != "0.0.0"
