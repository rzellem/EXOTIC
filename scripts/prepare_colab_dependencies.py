"""Prepare local Colab dependency wheels; never upload or publish them.

Only development files are removed. Runtime code, compiled extensions,
catalogues and licence files are retained. The wheel build tag distinguishes
these distributions without changing their upstream package versions.
"""
import argparse
from collections import Counter
from email.parser import BytesParser
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
from zipfile import ZipFile

from packaging.utils import parse_wheel_filename
from packaging.version import Version
from wheel.wheelfile import WheelFile


DEPENDENCIES = ('astroquery==0.4.11', 'photutils==3.1.1',
                'colour-science==0.4.7', 'colour-demosaicing==0.2.7')
SLIM_PACKAGES = {'astroquery', 'photutils', 'colour-science', 'colour-demosaicing'}
BUILD_TAG = '1colabslim'


def exclusion_reason(name):
    parts = PurePosixPath(name).parts
    if not parts or parts[0] not in {'astroquery', 'photutils', 'colour', 'colour_demosaicing'}:
        return None
    # Preserve all upstream notices, including notices inside test directories.
    if any(part.lower().startswith(('license', 'licence', 'copying', 'notice')) for part in parts):
        return None
    if 'tests' in parts:
        return 'tests and test fixtures'
    if parts[:2] == ('colour', 'htmlcov'):
        return 'HTML coverage reports'
    if len(parts) > 1 and parts[1] == 'examples':
        return 'examples'
    return None


def slim_wheel(source, output):
    source, output = Path(source), Path(output)
    package, version, _build, _tags = parse_wheel_filename(source.name)
    if package not in SLIM_PACKAGES:
        raise ValueError(f'No development-file rules defined for {package}')
    pieces = source.name.split('-')
    # The original distribution/version and platform tags remain unchanged.
    filename = '-'.join([pieces[0], pieces[1], BUILD_TAG, *pieces[-3:]])
    destination = output / filename
    removed = []
    retained = []
    with ZipFile(source) as archive, WheelFile(destination, 'w') as result:
        for member in archive.infolist():
            if member.is_dir() or member.filename.endswith('.dist-info/RECORD'):
                continue
            reason = exclusion_reason(member.filename)
            if reason:
                removed.append(dict(path=member.filename, reason=reason,
                                    compressed_bytes=member.compress_size, bytes=member.file_size))
                continue
            payload = archive.read(member)
            if member.filename.endswith('.dist-info/WHEEL'):
                # Match the build tag in the filename; other metadata is intact.
                payload = b'\n'.join(line for line in payload.split(b'\n')
                                     if not line.startswith(b'Build:'))
                payload = payload.rstrip(b'\n') + b'\nBuild: ' + BUILD_TAG.encode() + b'\n'
            result.writestr(member.filename, payload, compress_type=member.compress_type)
            retained.append(member.filename)
    # This software-distribution verification checks that retained files and
    # licences are unchanged. WheelFile regenerates RECORD once at this write
    # boundary for pip's required distribution integrity; no science is hashed.
    with ZipFile(source) as original, ZipFile(destination) as prepared:
        for name in retained:
            if not name.endswith('.dist-info/WHEEL'):
                if original.read(name) != prepared.read(name):
                    raise ValueError(f'Retained software file changed: {name}')
        metadata_name = next(n for n in retained if n.endswith('.dist-info/METADATA'))
        metadata = BytesParser().parsebytes(prepared.read(metadata_name))
        assert metadata['Version'] == str(version)
    return dict(package=str(package), version=str(version), source=source.name,
                wheel=filename, source_bytes=source.stat().st_size,
                wheel_bytes=destination.stat().st_size, retained_files=len(retained),
                removed_bytes=sum(item['bytes'] for item in removed),
                removed_files=removed,
                removed_file_counts=dict(Counter(item['reason'] for item in removed)),
                integrity='Retained payloads compared to upstream; RECORD regenerated for software installation integrity.')


def prepare_dependencies(output, upstream=None):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if upstream is None:
        upstream = output / 'upstream'
        upstream.mkdir(exist_ok=True)
        subprocess.check_call([
            sys.executable, '-m', 'pip', 'download', '--only-binary=:all:', '--no-deps',
            '--platform', 'manylinux_2_28_x86_64', '--platform', 'manylinux2014_x86_64',
            '--python-version', '3.13', '--implementation', 'cp', '--abi', 'cp313', '--abi', 'abi3',
            '--dest', str(upstream), *DEPENDENCIES])
    records = []
    selected = []
    for requirement in DEPENDENCIES:
        package, version = requirement.split('==')
        source = next(path for path in Path(upstream).glob('*.whl')
                      if parse_wheel_filename(path.name)[:2] == (package, Version(version)))
        record = slim_wheel(source, output)
        records.append(record)
        selected.append(output / record['wheel'])
    subprocess.check_call([sys.executable, '-m', 'pip', 'wheel', '--no-deps',
                           '--wheel-dir', str(output), 'imreg-dft==2.0.0'])
    imreg = next(output.glob('imreg_dft-2.0.0-*.whl'))
    with ZipFile(imreg) as archive:
        if not imreg.name.endswith('-none-any.whl') or any(n.endswith(('.so', '.pyd', '.dll')) for n in archive.namelist()):
            raise ValueError('Expected the pure Python imreg-dft wheel')
    selected.append(imreg)
    (output / 'dependency_audit.json').write_text(json.dumps(dict(
        slimmed=records, imreg_wheel=imreg.name,
        imreg_build='Pure Python upstream wheel; built once during distribution preparation.'), indent=2) + '\n', encoding='utf-8')
    return selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--upstream', type=Path, help='Use already downloaded upstream wheels')
    args = parser.parse_args()
    for wheel in prepare_dependencies(args.output, args.upstream):
        print(wheel)


if __name__ == '__main__':
    main()
