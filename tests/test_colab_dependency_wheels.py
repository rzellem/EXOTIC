"""Validate distribution trimming without changing runtime software files."""
from email.parser import BytesParser
import json
from pathlib import Path
import runpy
from zipfile import ZipFile

from packaging.utils import parse_wheel_filename
import pytest
from wheel.wheelfile import WheelFile

ROOT = Path(__file__).resolve().parents[1]
PREPARER = runpy.run_path(str(ROOT / 'scripts' / 'prepare_colab_dependencies.py'))


def test_slimming_retains_runtime_data_code_and_licences(tmp_path):
    upstream = tmp_path / 'astroquery-0.4.11-py3-none-any.whl'
    files = {
        'astroquery/__init__.py': b'__version__ = "0.4.11"\n',
        'astroquery/client.py': b'def query(): return 42\n',
        'astroquery/data/catalogue.fits': b'Necessary runtime data',
        'astroquery/tests/data/unused.fits': b'Test fixture',
        'astroquery/tests/test_client.py': b'assert True\n',
        'astroquery/tests/LICENSE': b'Fixture licence retained',
        'astroquery-0.4.11.dist-info/METADATA': b'Metadata-Version: 2.4\nName: astroquery\nVersion: 0.4.11\n',
        'astroquery-0.4.11.dist-info/WHEEL': b'Wheel-Version: 1.0\nRoot-Is-Purelib: true\nTag: py3-none-any\n',
        'astroquery-0.4.11.dist-info/licenses/LICENSE': b'Original licence',
    }
    with WheelFile(upstream, 'w') as archive:
        for name, payload in files.items():
            archive.writestr(name, payload)
    output = tmp_path / 'prepared'
    output.mkdir()
    report = PREPARER['slim_wheel'](upstream, output)
    assert report['removed_file_counts'] == {'tests and test fixtures': 2}
    prepared = output / report['wheel']
    assert parse_wheel_filename(prepared.name)[2] == (1, 'colabslim')
    with WheelFile(prepared) as archive:
        for name, payload in files.items():
            if name in ('astroquery/tests/data/unused.fits', 'astroquery/tests/test_client.py'):
                assert name not in archive.namelist()
            elif name.endswith('/WHEEL'):
                assert b'Build: 1colabslim' in archive.read(name)
            else:
                # Reading with WheelFile also checks its regenerated RECORD.
                assert archive.read(name) == payload


@pytest.mark.parametrize('name', ['colour/htmlcov/index.html', 'colour/examples/demo.py',
                                 'colour_demosaicing/examples/demo.ipynb',
                                 'photutils/psf/tests/data/image.fits'])
def test_development_file_rules(name):
    assert PREPARER['exclusion_reason'](name)


@pytest.mark.parametrize('name', ['colour/quality/datasets/spectrum.csv.gz',
                                 'colour_demosaicing/bayer/demosaicing/menon2007.py',
                                 'photutils/aperture/_batch_photometry.abi3.so',
                                 'barycorrpy/data/hip2.dat', 'astroquery/tests/data/LICENSE.txt'])
def test_runtime_files_and_notices_are_never_excluded(name):
    assert PREPARER['exclusion_reason'](name) is None


def test_prepared_dependencies_match_manifest_and_audit():
    directory = ROOT / 'wheelhouse' / 'colab'
    manifest = json.loads((directory / 'index.json').read_text())
    audit = json.loads((directory / 'dependency_audit.json').read_text())
    assert len(audit['slimmed']) == 4
    for record in audit['slimmed']:
        assert record['wheel'] in manifest['wheels']
        assert record['wheel_bytes'] < record['source_bytes']
        with WheelFile(directory / record['wheel']) as archive:
            assert not any(PREPARER['exclusion_reason'](name) for name in archive.namelist())
            metadata = BytesParser().parsebytes(archive.read(next(
                name for name in archive.namelist() if name.endswith('.dist-info/METADATA'))))
            assert metadata['Version'] == record['version']
            for name in archive.namelist():
                archive.read(name)
    assert audit['imreg_wheel'] in manifest['wheels']
    assert len(audit['remaining']) == 14
    assert all(item['wheel'] in manifest['wheels'] for item in audit['remaining'])
    with ZipFile(directory / audit['imreg_wheel']) as archive:
        assert 'imreg_dft/imreg.py' in archive.namelist()
        assert not any(name.endswith(('.so', '.dll', '.pyd')) for name in archive.namelist())
