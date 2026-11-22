"""Build local Linux EXOTIC/UltraNest wheels for reuse in Colab runtimes.

Run in a Linux build environment with build, auditwheel and patchelf installed.
This creates local artifacts; it does not publish or upload them.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

from packaging.utils import parse_wheel_filename


def write_colab_manifest(wheelhouse, wheels):
    """Expose this build's exact wheel versions to pip's standard selector."""
    wheels = sorted(wheels, key=lambda wheel: wheel.name)
    requirements = sorted({f'{name}=={version}' for wheel in wheels
                           for name, version, _build, _tags in [parse_wheel_filename(wheel.name)]})
    (wheelhouse / 'requirements.txt').write_text(
        '# Tested Colab software versions; native packages must use binary wheels.\n'
        + '\n'.join(requirements) + '\n', encoding='utf-8')
    (wheelhouse / 'index.json').write_text(
        json.dumps(dict(requirements=requirements, wheels=[wheel.name for wheel in wheels]), indent=2)
        + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path,
                        default=Path(__file__).resolve().parents[1],
                        help='EXOTIC checkout or source-distribution archive')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ultranest-version', help='Select a tested UltraNest version; default: latest available')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    workspace = output / datetime.now(timezone.utc).strftime('build-%Y%m%dT%H%M%S%fZ')
    raw = workspace / 'raw'
    raw.mkdir(parents=True)
    wheelhouse = output / 'wheelhouse'
    wheelhouse.mkdir(exist_ok=True)
    log_path = workspace / 'build.log'

    def run(command):
        command = list(map(str, command))
        print('RUN:', ' '.join(command), flush=True)
        completed = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        with log_path.open('a', encoding='utf-8') as log:
            log.write('RUN: ' + ' '.join(command) + '\n' + completed.stdout + '\n')
        print(completed.stdout, flush=True)
        completed.check_returncode()
        return completed.stdout

    run([sys.executable, '-m', 'pip', 'wheel', '--no-deps', '--wheel-dir', raw,
         args.source.resolve()])
    ultranest = 'ultranest' + ('==' + args.ultranest_version if args.ultranest_version else '')
    run([sys.executable, '-m', 'pip', 'wheel', '--no-deps', '--wheel-dir', raw, ultranest])
    run([sys.executable, '-m', 'pip', 'wheel', '--no-deps', '--wheel-dir', raw,
         'git+https://github.com/rzellem/exotic-colab.git'])
    records = []
    selected_wheels = []
    for wheel in sorted(raw.glob('*.whl')):
        if wheel.name.endswith('-none-any.whl'):
            shutil.copy2(wheel, wheelhouse / wheel.name)
            selected_wheels.append(wheelhouse / wheel.name)
            records.append(dict(wheel=wheel.name, compatible_tag='any',
                                analysis='Pure Python notebook display helper.'))
            continue
        analysis = run([sys.executable, '-m', 'auditwheel', 'show', wheel])
        match = re.search(r'platform tag:\s*"([^"\n]+)"', analysis)
        if match is None:
            raise RuntimeError(f'Could not read auditwheel platform analysis for {wheel.name}; see {log_path}')
        compatible_tag = match.group(1)
        repaired = workspace / 'repaired' / wheel.stem
        repaired.mkdir(parents=True)
        run([sys.executable, '-m', 'auditwheel', 'repair', '--plat', compatible_tag,
             '--wheel-dir', repaired, wheel])
        for built in sorted(repaired.glob('*.whl')):
            shutil.copy2(built, wheelhouse / built.name)
            selected_wheels.append(wheelhouse / built.name)
        records.append(dict(wheel=wheel.name, compatible_tag=compatible_tag,
                            analysis=analysis))
    from prepare_colab_dependencies import prepare_dependencies
    dependency_directory = workspace / 'dependencies'
    for prepared in prepare_dependencies(dependency_directory):
        shutil.copy2(prepared, wheelhouse / prepared.name)
        selected_wheels.append(wheelhouse / prepared.name)
    shutil.copy2(dependency_directory / 'dependency_audit.json', wheelhouse / 'dependency_audit.json')
    write_colab_manifest(wheelhouse, selected_wheels)
    report = dict(python=sys.version, platform=platform.platform(), source=str(args.source.resolve()),
                  build_workspace=str(workspace), wheels=records)
    (workspace / 'wheel_audit.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('Binary wheels ready:', wheelhouse, flush=True)


if __name__ == '__main__':
    main()
