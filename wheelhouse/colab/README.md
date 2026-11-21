# Prebuilt Colab software

These are the EXOTIC, UltraNest and notebook-display wheels validated in the
Colab CPU runtime on 22 November 2026. The complete set supports CPython 3.13,
Linux x86_64 and glibc 2.38 or newer. EXOTIC's ctypes library uses native ABI 6
and has no Python C ABI dependency; UltraNest's extension uses CPython 3.13.
The package versions and wheel filenames are recorded in `requirements.txt`
and `index.json`. The original wheel license files are preserved.

The `.numpy24` EXOTIC wheel revises the declared NumPy minimum to 2.4.6.
Its compiled library and native sources are identical to the previous tested
wheel; only dependency metadata and the matching version label changed. The
installer constrains NumPy to the runtime's existing installed version, so a
normal Colab setup keeps its preloaded NumPy instead of upgrading it.

The `.colabdeps` revision also accepts Colab's requests 2.32.4, Numba 0.64.0
and Bokeh 3.8.2. The installer reads this wheel's actual dependency requirements
and constrains compatible installed NumPy/requests/Numba/Bokeh versions, avoiding
the google-colab, PyTensor and Panel conflicts caused by unnecessary upgrades.
Missing or incompatible versions remain available for pip to resolve. Bootstrap
packages are installed together without upgrading already-satisfied versions.
Manifest retrieval, wheel downloads and pip installation are timed separately.
This revision changes only EXOTIC's dependency metadata and version label; all
64 other package files, including its native library, are unchanged.

The `.localmath` revision includes EXOTIC's local Python reference calculations
for orbital geometry, transit/eclipse integration and eclipse timing. Native
fits, Python comparisons and joint fits no longer import PyLightcurve or load
its catalogues. The original PyLightcurve 4.0.4 quadrature tables and numerical
rules are retained under its included MIT license. Only the optional
`exotethys` extra retains that external dependency for examples that call it.
Normal LDTK/PHOENIX stellar profiles are unchanged.

This revision repackages the tested Linux wheel with the revised `elca.py` and
`joint_fitter.py`, two local reference modules, dependency metadata and matching
version label. The compiled library, native sources and the other 60 package
files are byte-identical to the previous wheel. Wheel RECORD entries are
regenerated for software-distribution integrity; no science data is hashed.

The installer downloads the manifest's wheels and uses pip's binary-only
selection for these three packages. Its JSON manifest works with GitHub's raw
file MIME types, without requiring an HTML package-index service.
Ordinary dependencies continue to resolve through the normal package index.
If a future Colab Python ABI or platform needs another UltraNest wheel, pip
reports the incompatibility instead of silently compiling it from source.
If the retained compatible packages conflict with another dependency, pip
reports that conflict normally.

From an existing checkout, the wheelhouse can already be used locally:

```bash
python scripts/install_colab.py --wheelhouse wheelhouse/colab
```

For shared automatic downloads, the default installer URL points to this
directory on the existing `rzellem/EXOTIC` branch
`michael_fitzgeralds_wonderful_branch_of_magic`. Commit and push the prepared
files to make them available at that location; neither installer nor builder
uploads or publishes anything.

The original native wheel passed 818 Colab tests, including replay of 10,325
saved HAT-P-26b/WASP-16b posterior draws. The `.localmath` revision passed 3,090
local tests and 957 tests in a fresh Colab CPU runtime (Python 3.13.18, glibc
2.39, NumPy 2.4.6). These include frozen PyLightcurve 4.0.4 vectors, the joint
transit/RV/ephemeris likelihood, Python/native model comparisons, orbital replay,
photometry and notebook imports. The missing-library load-error fixture is
isolated from the source builder so it does not deliberately compile a library.

The fresh runtime downloaded all three wheels automatically from this branch's
published JSON manifest, retained its original NumPy, and installed/imported no
PyLightcurve. Installation took 91.9 seconds and the notebook imports took 5.1
seconds. Normal native loads and the regression suite required no C++ compiler
calls or ZIP uploads. Drive mounting and full science-image reductions were not
exercised in this setup test. Pip reported existing dependency conflicts with
the runtime's google-colab/requests, pytensor/numba and panel/bokeh requirements;
the EXOTIC checks above passed despite those warnings.

Rebuild the wheel set with `scripts/build_colab_wheels.py`. It also regenerates
the exact version requirements and wheel manifest from that build's outputs.
Do not reuse an old EXOTIC wheel after changing its Python or native source;
the complete EXOTIC wheel is a snapshot of both.
