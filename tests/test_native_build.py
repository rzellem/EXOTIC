"""Automatic compilation, reuse, failures and contention between bulk workers."""
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import platform
import runpy
import subprocess
import sys
import time

import pytest

from exotic.api import native_build, limbdarkfast_cpp


def fake_compiler(monkeypatch):
    calls = []

    def compile(source, library, build):
        calls.append(library)
        time.sleep(0.03)
        library.write_bytes(b'native library')
        return dict(abi_version=native_build.ABI_VERSION, sys_platform=sys.platform,
                    machine=platform.machine())

    monkeypatch.setattr(native_build, '_compile', compile)
    return calls


def test_missing_build_is_created_once_and_reused(tmp_path, monkeypatch):
    calls = fake_compiler(monkeypatch)
    first = native_build.ensure_hybrid4_native(output_dir=tmp_path)
    assert native_build.ensure_hybrid4_native(output_dir=tmp_path) == first
    assert calls == [first]


def test_concurrent_bulk_workers_compile_once(tmp_path, monkeypatch):
    calls = fake_compiler(monkeypatch)
    with ThreadPoolExecutor(max_workers=10) as pool:
        libraries = list(pool.map(lambda _: native_build.ensure_hybrid4_native(output_dir=tmp_path), range(10)))
    assert len(calls) == 1
    assert len(set(libraries)) == 1


def test_previous_native_abi_is_rebuilt_automatically(tmp_path, monkeypatch):
    calls = fake_compiler(monkeypatch)
    native_build.ensure_hybrid4_native(output_dir=tmp_path)
    metadata_path = tmp_path / 'build_info.json'
    metadata = json.loads(metadata_path.read_text())
    metadata['abi_version'] = native_build.ABI_VERSION - 1
    metadata_path.write_text(json.dumps(metadata))
    native_build.ensure_hybrid4_native(output_dir=tmp_path)
    assert len(calls) == 2
    assert json.loads(metadata_path.read_text())['abi_version'] == native_build.ABI_VERSION


def test_wheel_extraction_times_do_not_require_a_compiler(tmp_path, monkeypatch):
    calls = fake_compiler(monkeypatch)
    library = native_build.ensure_hybrid4_native(output_dir=tmp_path)
    metadata_path = tmp_path / 'build_info.json'
    metadata = json.loads(metadata_path.read_text())
    metadata['packaged'] = True
    metadata_path.write_text(json.dumps(metadata))
    # A wheel can be extracted with its header newer than the compiled library.
    header = tmp_path / 'limbdarkfast.hpp'
    header.write_text('installed header')
    stamp = library.stat().st_mtime_ns + 1_000_000_000
    os.utime(header, ns=(stamp, stamp))
    native_build.ensure_hybrid4_native(source_dir=tmp_path, output_dir=tmp_path)
    assert len(calls) == 1


def test_changed_sources_rebuild_and_force_rebuilds(tmp_path, monkeypatch):
    source = tmp_path / 'source'
    source.mkdir()
    for name in ('hybrid4_capi.cpp', 'limbdarkfast.hpp', 'claret.hpp', 'photometry.hpp', 'likelihood.hpp', 'image_preparation.hpp'):
        (source / name).write_text('source')
    calls = fake_compiler(monkeypatch)
    library = native_build.ensure_hybrid4_native(source_dir=source, output_dir=tmp_path)
    stamp = library.stat().st_mtime_ns + 1_000_000_000
    os.utime(source / 'limbdarkfast.hpp', ns=(stamp, stamp))
    native_build.ensure_hybrid4_native(source_dir=source, output_dir=tmp_path)
    native_build.ensure_hybrid4_native(source_dir=source, output_dir=tmp_path, force=True)
    assert len(calls) == 3


def test_compiler_failure_propagates_without_accepting_partial_binary(tmp_path, monkeypatch):
    def fail(source, library, build):
        library.write_bytes(b'incomplete')
        raise subprocess.CalledProcessError(2, ['compiler'])
    monkeypatch.setattr(native_build, '_compile', fail)
    with pytest.raises(subprocess.CalledProcessError):
        native_build.ensure_hybrid4_native(output_dir=tmp_path)
    assert not native_build._is_current(tmp_path / native_build.library_filename(), native_build.NATIVE)
    calls = fake_compiler(monkeypatch)
    native_build.ensure_hybrid4_native(output_dir=tmp_path)
    assert len(calls) == 1


def test_cpp_first_use_prepares_backend_and_reports_build_failure(monkeypatch):
    limbdarkfast_cpp._library.cache_clear()
    def fail(**kwargs):
        raise RuntimeError('C++ compiler is missing')
    monkeypatch.setattr(limbdarkfast_cpp, 'ensure_hybrid4_native', fail)
    with pytest.raises(RuntimeError, match=r'C\+\+ compiler is missing'):
        limbdarkfast_cpp.hybrid4_flux([0.1], 0.1, [0.3, 0.2, 0.1, 0.05])
    limbdarkfast_cpp._library.cache_clear()


@pytest.mark.parametrize('inplace', [False, True])
def test_package_build_supplies_binary_and_metadata(tmp_path, monkeypatch, inplace):
    from setuptools import Distribution, Extension
    module = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'setup.py'))
    builder = module['native_build']
    source_output = tmp_path / 'source/exotic/api/native'
    source_output.mkdir(parents=True)
    monkeypatch.setattr(builder, 'NATIVE', source_output)
    calls = fake_compiler(monkeypatch)
    monkeypatch.setattr(builder, '_compile', native_build._compile)
    dist = Distribution(dict(packages=['exotic', 'exotic.api'],
                             package_dir={'': str(tmp_path / 'source')},
                             ext_modules=[Extension('exotic.api.native.limbdarkfast', sources=[])]))
    cmd = module['Hybrid4BuildExt'](dist)
    cmd.build_lib = str(tmp_path / 'built')
    cmd.build_temp = str(tmp_path / 'obj')
    cmd.inplace = inplace
    cmd.ensure_finalized()
    cmd.run()
    built = tmp_path / 'built/exotic/api/native' / builder.library_filename()
    assert built.is_file()
    assert str(built.parent / 'build_info.json') in cmd.get_outputs()
    if inplace:
        assert (source_output / builder.library_filename()).is_file()
        assert json.loads((source_output / 'build_info.json').read_text())['abi_version'] == native_build.ABI_VERSION
    assert len(calls) == 1


@pytest.mark.parametrize('platform_tag', ['linux_x86_64', 'win_amd64', 'macosx_11_0_arm64'])
def test_ctypes_wheel_keeps_platform_without_a_python_abi(platform_tag, monkeypatch):
    from setuptools import Distribution, Extension
    from setuptools.command.bdist_wheel import bdist_wheel
    module = runpy.run_path(str(Path(__file__).resolve().parents[1] / 'setup.py'))
    monkeypatch.setattr(bdist_wheel, 'get_tag', lambda self: ('cp313', 'cp313', platform_tag))
    command = module['Hybrid4BdistWheel'](Distribution(dict(
        ext_modules=[Extension('exotic.api.native.limbdarkfast', sources=[])])))
    assert command.get_tag() == ('py3', 'none', platform_tag)
    assert command.distribution.has_ext_modules()
