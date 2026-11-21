"""Resolve the real Colab wheels through pip without installing into this host."""
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from email.parser import BytesParser
import json
from pathlib import Path
import runpy
import subprocess
from threading import Thread
from urllib.error import HTTPError
from zipfile import ZipFile

from packaging.requirements import Requirement

import pytest


ROOT = Path(__file__).resolve().parents[1]
WHEELHOUSE = ROOT / 'wheelhouse' / 'colab'


@pytest.fixture
def wheel_server():
    requests = []

    class QuietHandler(SimpleHTTPRequestHandler):
        # GitHub raw files are not served as HTML package indexes. The JSON
        # manifest and wheel downloads must work with plain MIME types.
        extensions_map = {'.json': 'text/plain', '.whl': 'application/octet-stream'}

        def do_GET(self):
            requests.append(self.path)
            super().do_GET()

        def log_message(self, *_args):
            pass

    server = ThreadingHTTPServer(('127.0.0.1', 0), partial(QuietHandler, directory=str(WHEELHOUSE)))
    worker = Thread(target=server.serve_forever)
    worker.start()
    try:
        yield dict(url=f'http://127.0.0.1:{server.server_port}', requests=requests)
    finally:
        server.shutdown()
        server.server_close()
        worker.join()


def resolve_with_pip(monkeypatch, tmp_path, source, python_version='3.13'):
    installer = runpy.run_path(str(ROOT / 'scripts' / 'install_colab.py'))
    report = tmp_path / 'pip-report.json'
    output = []

    def dry_run(command):
        # Execute the actual install command, adding pip's foreign-runtime
        # resolver options so Linux wheels can be tested on this Windows host.
        completed = subprocess.run([
            *command, '--dry-run', '--ignore-installed', '--no-deps', '--no-index',
            '--platform', 'manylinux_2_38_x86_64', '--platform', 'manylinux_2_17_x86_64',
            '--platform', 'manylinux2014_x86_64', '--implementation', 'cp',
            '--python-version', python_version, '--abi', 'cp' + python_version.replace('.', ''),
            '--report', str(report),
        ], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output.append(completed.stdout)
        completed.check_returncode()

    monkeypatch.setattr(installer['subprocess'], 'check_call', dry_run)
    installer['install_prebuilt'](source)
    return json.loads(report.read_text(encoding='utf-8')), ''.join(output)


@pytest.mark.parametrize('remote', [False, True], ids=['local-directory', 'automatic-download'])
def test_real_colab_wheels_resolve_without_source_builds(monkeypatch, tmp_path, wheel_server, remote):
    source = wheel_server['url'] if remote else WHEELHOUSE
    report, output = resolve_with_pip(monkeypatch, tmp_path, source)
    installed = {item['metadata']['name'].lower(): item for item in report['install']}
    assert set(installed) == {'exotic', 'ultranest', 'exoticcolab'}
    assert installed['ultranest']['metadata']['version'] == '4.6.3'
    assert all(item['download_info']['url'].endswith('.whl') for item in installed.values())
    assert 'Building wheel' not in output
    if remote:
        assert '/index.json' in wheel_server['requests']
        assert all('/' + wheel.name in wheel_server['requests']
                   for wheel in WHEELHOUSE.glob('*.whl'))


def test_incompatible_python_does_not_start_a_source_build(monkeypatch, tmp_path, wheel_server):
    # UltraNest's CPython 3.13 extension cannot be loaded by CPython 3.14.
    with pytest.raises(subprocess.CalledProcessError) as failure:
        resolve_with_pip(monkeypatch, tmp_path, wheel_server['url'], python_version='3.14')
    assert 'No matching distribution found for ultranest==4.6.3' in failure.value.output
    assert 'Building wheel' not in failure.value.output


def test_missing_download_is_reported_without_a_source_retry(monkeypatch, tmp_path, wheel_server):
    with pytest.raises(HTTPError, match='404'):
        resolve_with_pip(monkeypatch, tmp_path, wheel_server['url'] + '/not-published')
    assert not (tmp_path / 'pip-report.json').exists()


def test_new_build_index_does_not_advertise_previous_builds(tmp_path):
    builder = runpy.run_path(str(ROOT / 'scripts' / 'build_colab_wheels.py'))
    current = list(WHEELHOUSE.glob('*.whl'))
    stale = tmp_path / 'exotic-99.0-py3-none-any.whl'
    stale.touch()
    builder['write_colab_manifest'](tmp_path, current)
    index = (tmp_path / 'index.json').read_text(encoding='utf-8')
    requirements = (tmp_path / 'requirements.txt').read_text(encoding='utf-8')
    assert stale.name not in index
    assert 'exotic==99.0' not in requirements
    assert all(wheel.name in index for wheel in current)
    assert 'ultranest==4.6.3' in requirements


def test_bundled_exotic_accepts_colabs_numpy_version():
    manifest = json.loads((WHEELHOUSE / 'index.json').read_text(encoding='utf-8'))
    wheel = next(name for name in manifest['wheels'] if name.startswith('exotic-'))
    with ZipFile(WHEELHOUSE / wheel) as archive:
        metadata_path = next(name for name in archive.namelist() if name.endswith('.dist-info/METADATA'))
        metadata = BytesParser().parsebytes(archive.read(metadata_path))
    numpy_requirement = next(Requirement(value) for value in metadata.get_all('Requires-Dist')
                             if Requirement(value).name == 'numpy')
    assert '2.4.6' in numpy_requirement.specifier
    assert metadata['Version'] in '\n'.join(manifest['requirements'])
