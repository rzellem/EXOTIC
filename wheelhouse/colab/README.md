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
Likewise, if another dependency cannot use the retained NumPy version, pip
reports that dependency conflict instead of replacing NumPy.

From an existing checkout, the wheelhouse can already be used locally:

```bash
python scripts/install_colab.py --wheelhouse wheelhouse/colab
```

For shared automatic downloads, the default installer URL points to this
directory on the existing `rzellem/EXOTIC` branch
`michael_fitzgeralds_wonderful_branch_of_magic`. Commit and push the prepared
files to make them available at that location; neither installer nor builder
uploads or publishes anything.

The actual Colab wheel installation passed 818 tests. Its normal native loads
and replay of 10,325 saved HAT-P-26b/WASP-16b posterior draws made zero compiler
calls. That validation was saved-data replay, not a new full image reduction.
The automatic-download path is tested locally through a loopback HTTP server
with pip's CPython 3.13 / Linux resolver; it has not been tested against a
published repository URL.

Rebuild the wheel set with `scripts/build_colab_wheels.py`. It also regenerates
the exact version requirements and wheel manifest from that build's outputs.
Do not reuse an old EXOTIC wheel after changing its Python or native source;
the complete EXOTIC wheel is a snapshot of both.
