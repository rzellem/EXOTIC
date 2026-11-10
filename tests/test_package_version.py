"""Package version detection for installed distributions and source checkouts."""

from importlib.util import module_from_spec, spec_from_file_location
from importlib import metadata
from pathlib import Path

from exotic.version import __version__ as source_version


PACKAGE_INIT = Path(__file__).resolve().parents[1] / 'exotic' / '__init__.py'


def load_package():
    spec = spec_from_file_location('exotic', PACKAGE_INIT)
    package = module_from_spec(spec)
    spec.loader.exec_module(package)
    return package


def test_installed_distribution_version_takes_precedence(monkeypatch):
    monkeypatch.setattr(metadata, 'version', lambda name: '9.8.7')

    package = load_package()

    assert package.__version__ == '9.8.7'


def test_source_checkout_uses_generated_version_without_distribution_metadata(monkeypatch):
    def not_installed(name):
        raise metadata.PackageNotFoundError(name)

    monkeypatch.setattr(metadata, 'version', not_installed)

    package = load_package()

    assert package.__version__ == source_version
