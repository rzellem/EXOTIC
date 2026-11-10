# Local Windows wheels

These wheels target standard Windows x64 CPython. UltraNest 4.5.0 is included
for Python 3.12, 3.13, and 3.14:

- `ultranest-4.5.0-cp312-cp312-win_amd64.whl` — CPython 3.12
- `ultranest-4.5.0-cp313-cp313-win_amd64.whl` — CPython 3.13
- `ultranest-4.5.0-cp314-cp314-win_amd64.whl` — CPython 3.14

Python 3.14 also needs a local SEP wheel for Astroalign's dependency:

- `sep-1.4.1-cp314-cp314-win_amd64.whl` — CPython 3.14

SEP 1.4.1 has no upstream Windows CPython 3.14 wheel as of 2026-10-24. Including
this build allows pip to install it without requiring users to compile SEP.
Python 3.12/3.13 can use SEP wheels from PyPI. Astropy 8 uses a stable-ABI wheel
that works on Python 3.14, so it does not need a separate local `cp314` wheel.

They were built locally from the packages' source distributions with the
Visual Studio 2022 C++ toolchain. The new CPython 3.14 builds used Python 3.14.9,
Cython 3.3.0, and NumPy 2.5.4. The compiled UltraNest extensions were imported
under their matching interpreters. The 3.14 environment also passed native
pair-distance and SEP source-extraction tests and a small UltraNest sampling run.

SHA-256:

These recorded hashes identify the compiled software distributed to Windows
users. Filenames and ABI tags describe compatibility but do not distinguish
altered wheel contents; they do not replace this software-distribution identity.

- CPython 3.12: `69C57FC84FD475AEA6A959F473F942192A6839860849D626456E0F61406B06FF`
- CPython 3.13: `09CA8F766E039F072B9182B6072F0E318FA6173EA516368B988E5F9F56A9834C`
- UltraNest CPython 3.14: `3A49E3BC905E55000625A791F12CB692CD4A3342E11A38A0C5F73287DC731644`
- SEP CPython 3.14: `C273D9664381306B76CDF4F85873337979812964C9E63FC7F1182CE98766ECD2`

The new 3.14 digests above reuse pip's recorded wheel-build results; the
unchanged accepted wheels are not repeatedly scanned.

`run_exotic_windows.bat` adds this directory as a pip `--find-links` source and
requires a binary UltraNest distribution, preventing a fallback source build.
It uses the active environment's Python for both `python -m pip` and the GUI.

## Rebuild the Python 3.14 wheels

In a Visual Studio **x64 Native Tools Command Prompt**, activate a Python 3.14
virtual environment, then run from this repository's root:

```bat
python -m pip install --upgrade pip setuptools wheel cython numpy
python -m pip wheel --no-deps --no-binary=ultranest ultranest==4.5.0 --wheel-dir wheelhouse
```

For SEP 1.4.1, apply the included [license-files patch](sep-1.4.1-license-files.patch)
to its extracted source before building. It adds the upstream LGPL, BSD, and MIT
license texts to the wheel through setuptools' packaging configuration:

```bat
mkdir tmp\sep-build
python -m pip download --no-deps --no-binary=sep sep==1.4.1 --dest tmp\sep-build
tar -xf tmp\sep-build\sep-1.4.1.tar.gz -C tmp\sep-build
git apply --directory=tmp/sep-build/sep-1.4.1 wheelhouse/sep-1.4.1-license-files.patch
python -m pip wheel --no-deps tmp\sep-build\sep-1.4.1 --wheel-dir wheelhouse
```

Original sources: [UltraNest 4.5.0](https://pypi.org/project/ultranest/4.5.0/#files)
and [SEP 1.4.1](https://pypi.org/project/sep/1.4.1/#files). The SEP patch changes
packaging metadata to include the original license files.

Pip's isolated build environments install each package's declared build
dependencies. Retest the compiled extensions after rebuilding, and record the
new hashes emitted by pip at that build boundary. The compatibility workflow
installs the matching bundled UltraNest wheel on each Windows Python version
and the bundled SEP wheel on Python 3.14, then runs the full EXOTIC test suite
with minimum and latest Astropy versions.
