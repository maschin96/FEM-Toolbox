"""Smoke tests for the installable package boundary."""

from importlib import import_module
from importlib.metadata import metadata, version


def test_package_can_be_imported() -> None:
    """The public package namespace is importable after installation."""
    package = import_module("fem_toolbox")

    assert package.__name__ == "fem_toolbox"


def test_distribution_exposes_release_version() -> None:
    """Installed metadata identifies the M11 development version."""
    distribution = metadata("fem-toolbox")

    assert distribution["Name"] == "fem-toolbox"
    assert version("fem-toolbox") == "1.6.0.dev1"
    assert import_module("fem_toolbox").__version__ == version("fem-toolbox")
