"""Setuptools entry point; package metadata lives in setup.cfg and pyproject.toml.

Build distributions with ``python -m build``. Versioning is generated from Git
by setuptools_scm; publishing is handled by the Python Package workflow.
"""
from importlib.util import module_from_spec, spec_from_file_location
import json
import os
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.bdist_wheel import bdist_wheel


ROOT = Path(__file__).resolve().parent
spec = spec_from_file_location('exotic_native_build', ROOT / 'exotic/api/native_build.py')
native_build = module_from_spec(spec)
spec.loader.exec_module(native_build)


class Hybrid4BuildExt(build_ext):
    """Ship the C ABI library as a platform-specific extension without PyInit."""

    def get_ext_filename(self, name):
        return os.path.join(*name.split('.')[:-1], native_build.library_filename())

    def build_extension(self, ext):
        destination = Path(self.get_ext_fullpath(ext.name)).resolve().parent
        native_build.ensure_hybrid4_native(
            output_dir=destination, build_dir=Path(self.build_temp) / 'hybrid4',
            force=bool(self.force),
        )
        metadata_path = destination / 'build_info.json'
        metadata = json.loads(metadata_path.read_text(encoding='utf-8'))
        metadata['packaged'] = True
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')

    def copy_extensions_to_source(self):
        super().copy_extensions_to_source()
        # Metadata travels with the binary so first use need not recompile it.
        for ext in self.extensions:
            built = Path(self.build_lib) / Path(self.get_ext_filename(ext.name))
            self.copy_file(str(built.parent / 'build_info.json'),
                           str(native_build.NATIVE / 'build_info.json'))

    def get_outputs(self):
        outputs = super().get_outputs()
        return outputs + [str(Path(path).parent / 'build_info.json') for path in outputs]


class Hybrid4BdistWheel(bdist_wheel):
    """The ctypes C library depends on the OS/CPU, not the Python C ABI."""

    def get_tag(self):
        _python, _abi, platform_tag = super().get_tag()
        return 'py3', 'none', platform_tag


if __name__ == "__main__":
    setup(
        ext_modules=[Extension('exotic.api.native.limbdarkfast',
                               sources=['exotic/api/native/hybrid4_capi.cpp'],
                               depends=['exotic/api/native/limbdarkfast.hpp', 'exotic/api/native/claret.hpp', 'exotic/api/native/photometry.hpp', 'exotic/api/native/likelihood.hpp', 'exotic/api/native/image_preparation.hpp'])],
        cmdclass={'build_ext': Hybrid4BuildExt, 'bdist_wheel': Hybrid4BdistWheel},
    )
