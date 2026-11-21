"""Build hybrid4 for installs and checkouts, using only the standard library.

Concurrent bulk workers share an OS lock. Source modification times determine
when to rebuild; no science files are read or hashed.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import argparse
import errno
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


NATIVE = Path(__file__).resolve().parent / 'native'
ABI_VERSION = 6


def library_filename():
    if sys.platform == 'win32':
        return 'limbdarkfast.dll'
    return 'liblimbdarkfast.dylib' if sys.platform == 'darwin' else 'liblimbdarkfast.so'


def msvc_environment():
    """Discover compiler paths without printing the compiler environment."""
    locator = shutil.which('vswhere') or r'C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe'
    if not Path(locator).is_file():
        raise RuntimeError('C++ compilation needs MSVC C++ Build Tools or an available cl.exe.')
    installation = subprocess.check_output([
        locator, '-latest', '-products', '*', '-requires',
        'Microsoft.VisualStudio.Component.VC.Tools.x86.x64', '-property', 'installationPath',
    ], text=True).strip()
    if not installation:
        raise RuntimeError('MSVC C++ Build Tools are not installed.')
    setup = Path(installation) / 'VC' / 'Auxiliary' / 'Build' / 'vcvars64.bat'
    result = subprocess.check_output(f'cmd.exe /d /s /c ""{setup}" >nul && set"', text=True)
    # Compiler discovery paths, not EXOTIC application configuration.
    compiler_environment = {name.upper(): value for name, value in os.environ.items()}
    for line in result.splitlines():
        name, separator, value = line.partition('=')
        if separator and name and not name.startswith('='):
            compiler_environment[name.upper()] = value
    return compiler_environment


@contextmanager
def _build_lock(path):
    # OS locks release on process exit; there is no stale-lock timeout.
    with path.open('a+b') as stream:
        if sys.platform == 'win32':
            import msvcrt
            if stream.seek(0, os.SEEK_END) == 0:
                stream.write(b'\0')
                stream.flush()
            while True:
                stream.seek(0)
                try:
                    msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError as exc:
                    if exc.errno not in (errno.EACCES, errno.EAGAIN, errno.EDEADLK):
                        raise
                    time.sleep(0.1)
            try:
                yield
            finally:
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            fcntl.flock(stream, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(stream, fcntl.LOCK_UN)


def _is_current(library, source_dir):
    try:
        stamp = library.stat()
        metadata = json.loads((library.parent / 'build_info.json').read_text(encoding='utf-8'))
        if not (stamp.st_size > 0 and metadata['abi_version'] == ABI_VERSION
                and metadata['sys_platform'] == sys.platform
                and metadata['machine'] == platform.machine()):
            return False
        # Wheel installers assign new extraction times. A packaged binary is
        # already built with its bundled sources; those times are not edits.
        # Editable/source checkouts still rebuild after real source changes.
        source_root = source_dir.parents[2]
        checkout = (source_root / 'setup.py').is_file() and (source_root / 'pyproject.toml').is_file()
        if metadata.get('packaged') and not checkout:
            return True
        return all(path.stat().st_mtime_ns <= stamp.st_mtime_ns for path in (
            source_dir / 'hybrid4_capi.cpp', source_dir / 'limbdarkfast.hpp',
            source_dir / 'claret.hpp', source_dir / 'photometry.hpp',
            source_dir / 'likelihood.hpp', source_dir / 'image_preparation.hpp', Path(__file__)))
    except (OSError, ValueError, KeyError):
        return False


def _compile(source_dir, library, build):
    source = source_dir / 'hybrid4_capi.cpp'
    compiler_environment = None
    if sys.platform == 'win32':
        compiler_environment = msvc_environment() if shutil.which('cl') is None else dict(os.environ)
        compiler = shutil.which('cl', path=compiler_environment.get('PATH', compiler_environment.get('Path')))
        if not compiler:
            raise RuntimeError('MSVC setup did not provide cl.exe.')
        command = [compiler, '/nologo', '/O2', '/fp:precise', '/EHsc', '/std:c++14', '/LD',
                   str(source), f'/Fo{build / "hybrid4_capi.obj"}', '/link',
                   f'/OUT:{library}', f'/IMPLIB:{build / "limbdarkfast.lib"}']
    else:
        compiler = shutil.which('c++') or shutil.which('clang++') or shutil.which('g++')
        if not compiler:
            raise RuntimeError('C++ compilation needs an installed c++, clang++ or g++ compiler.')
        command = [compiler, '-std=c++11', '-O3', '-fPIC',
                   '-dynamiclib' if sys.platform == 'darwin' else '-shared', str(source), '-o', str(library)]
    print('Compiling hybrid4 with', compiler, flush=True)
    subprocess.run(command, cwd=build, env=compiler_environment, check=True)
    return {
        'compiler': compiler, 'command': command, 'platform': platform.platform(),
        'sys_platform': sys.platform, 'machine': platform.machine(), 'abi_version': ABI_VERSION,
        'built_utc': datetime.now(timezone.utc).isoformat(),
        'upstream_revision': '120caecc7d6d1d60c10cf0ed4ec27e4dede952db',
        'floating_point': 'Double precision; fast-math disabled', 'library': str(library),
    }


def ensure_hybrid4_native(*, source_dir=None, output_dir=None, build_dir=None, force=False):
    """Build when absent or stale; return the library path. Errors propagate."""
    source_dir = Path(source_dir or NATIVE).resolve()
    output_dir = Path(output_dir or NATIVE).resolve()
    library = output_dir / library_filename()
    if not force and _is_current(library, source_dir):
        return library
    build = Path(build_dir).resolve() if build_dir is not None else output_dir / 'build'
    output_dir.mkdir(parents=True, exist_ok=True)
    build.mkdir(parents=True, exist_ok=True)
    with _build_lock(build / 'hybrid4.lock'):
        if not force and _is_current(library, source_dir):
            return library
        (output_dir / 'build_info.json').unlink(missing_ok=True)
        metadata = _compile(source_dir, library, build)
        (output_dir / 'build_info.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
        print('Built', library, flush=True)
    return library


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--force', action='store_true', help='Rebuild even if current')
    args = parser.parse_args()
    print('hybrid4 ready:', ensure_hybrid4_native(force=args.force), flush=True)


if __name__ == '__main__':
    main()
