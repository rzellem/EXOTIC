"""Install the tested Colab wheels without compiling EXOTIC or UltraNest.

The default URL becomes available when the prepared repository files are
published. This installer only reads that location; it never uploads files.
Use --wheelhouse PATH to install from an existing local wheelhouse instead.
"""
import argparse
from email.parser import BytesParser
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import subprocess
import sys
import time
from tempfile import TemporaryDirectory
from urllib.request import urlopen, urlretrieve
from zipfile import ZipFile

from packaging.requirements import Requirement


WHEELHOUSE_URL = (
    'https://raw.githubusercontent.com/rzellem/EXOTIC/'
    'michael_fitzgeralds_wonderful_branch_of_magic/wheelhouse/colab'
)


def compatible_runtime_versions(exotic_wheel):
    """Keep Colab's shared packages when they satisfy this EXOTIC wheel."""
    with ZipFile(exotic_wheel) as archive:
        metadata = BytesParser().parsebytes(archive.read(next(
            name for name in archive.namelist() if name.endswith('.dist-info/METADATA'))))
    retained = {}
    for value in metadata.get_all('Requires-Dist', []):
        requirement = Requirement(value)
        if requirement.name.lower() not in ('numpy', 'requests', 'numba', 'bokeh'):
            continue
        if requirement.marker and not requirement.marker.evaluate({'extra': ''}):
            continue
        try:
            installed = version(requirement.name)
        except PackageNotFoundError:
            continue
        if installed in requirement.specifier:
            retained[requirement.name] = installed
            print(f'Keeping the runtime {requirement.name} {installed}.', flush=True)
        else:
            print(f'Runtime {requirement.name} {installed} does not satisfy {requirement}; '
                  'pip will resolve a compatible version.', flush=True)
    return retained


def install_prebuilt(wheelhouse=WHEELHOUSE_URL):
    started = time.perf_counter()
    timings = dict(manifest_seconds=0., wheel_downloads=[])
    source = str(wheelhouse).rstrip('/')
    with TemporaryDirectory(prefix='exotic-colab-') as temporary:
        if source.startswith(('https://', 'http://')):
            stage_started = time.perf_counter()
            with urlopen(source + '/index.json') as response:
                manifest = json.load(response)
            timings['manifest_seconds'] = time.perf_counter() - stage_started
            directory = Path(temporary)
            for wheel in manifest['wheels']:
                stage_started = time.perf_counter()
                urlretrieve(source + '/' + wheel, directory / wheel)
                timings['wheel_downloads'].append(dict(
                    wheel=wheel, seconds=time.perf_counter() - stage_started,
                    bytes=(directory / wheel).stat().st_size))
        else:
            directory = Path(source).resolve()
            manifest = json.loads((directory / 'index.json').read_text(encoding='utf-8'))
        exotic_wheel = directory / next(name for name in manifest['wheels'] if name.startswith('exotic-'))
        retained = compatible_runtime_versions(exotic_wheel)
        timings['retained_runtime_packages'] = retained
        constraints = []
        if retained:
            constraint_path = Path(temporary) / 'runtime-packages.txt'
            constraint_path.write_text(''.join(f'{name}=={value}\n' for name, value in retained.items()),
                                       encoding='utf-8')
            constraints = ['--constraint', str(constraint_path)]
        stage_started = time.perf_counter()
        subprocess.check_call([
            sys.executable, '-m', 'pip', 'install', '--upgrade',
            '--only-binary=exotic,ultranest,exoticcolab',
            '--find-links', str(directory), *constraints, *manifest['requirements'],
        ])
        timings['pip_seconds'] = time.perf_counter() - stage_started
    timings['total_seconds'] = time.perf_counter() - started
    print('Prebuilt EXOTIC, UltraNest and notebook display helper installed.', flush=True)
    print('Colab install timings: ' + json.dumps(timings), flush=True)
    return timings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheelhouse', default=WHEELHOUSE_URL,
                        help='Wheelhouse directory or URL; default: the EXOTIC Colab wheelhouse')
    args = parser.parse_args()
    install_prebuilt(args.wheelhouse)


if __name__ == '__main__':
    main()
