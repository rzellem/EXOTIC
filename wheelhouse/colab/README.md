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

The `.scipy116` revision accepts and retains Colab's SciPy 1.16.3 as well.
Its only wheel changes are the SciPy dependency minimum and matching version
label. All 64 other package files and the compiled native backend are unchanged.
Other environments may continue using newer compatible SciPy versions.

The fresh `.scipy116` runtime retained SciPy 1.16.3 alongside NumPy 2.4.6 and
the other Colab dependency versions. All 962 installed-wheel checks passed and
`pip check` reported no broken requirements. The full local suite under an
isolated SciPy 1.16.3 installation passed 3,095 tests, with the optional
Colab/IPython image check skipped locally and exercised successfully in Colab.

Normal setup took 55.7 seconds versus 73.2 seconds in the previous fresh runtime
(17.6 seconds, or 24%, shorter in this pair). Dependency payloads fell from
116.9 MB to 81.6 MB by retaining the existing 35.3 MB SciPy wheel installation.
Pip installation fell from 50.3 to 37.9 seconds; its unpack/replace/install phase
fell from 34.1 to 23.2 seconds. Notebook imports took 7.6 seconds. These are
individual fresh-runtime observations, with the same timing scope as above.

Scientific validation refitted the recorded HAT-P-26b (241 points) and WASP-16b
(72 points) Claret and hybrid4 inputs under both SciPy 1.16.3 and 1.18.1, using
three matched seeds: 24 independent fits and 12 version pairs. Photometry,
atmosphere coefficients, exposure integration, masks, priors and normal sampling
controls were held fixed. All pairs had identical models, posterior medians,
68% intervals, reported errors, chi-square, log evidence, prior bounds and
likelihood-call counts. Replaying 257 prior/likelihood points per input case
also gave zero numerical differences. The SciPy 1.18.1 reference was installed
into a separate validation directory and did not replace the runtime's SciPy;
that reference installation and all refits are excluded from setup timing.
This check reused saved photometry and did not rerun image reductions or
regenerate stellar atmosphere profiles. Native C++ compilation and ZIP uploads
were again unnecessary.

The `.colabstack` revision also accepts and retains Colab's Astropy 7.2.3,
Matplotlib 3.10.0 and scikit-image 0.25.2. Only these three dependency minima
and the matching wheel version label change; the other 64 package files,
including the compiled native backend, are unchanged. The installer pins all
eight shared packages only when the installed versions satisfy the wheel's
requirements. Newer compatible versions remain supported. Fresh-runtime
timings and matched fit comparisons are pending this revision's manual push.

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
the EXOTIC checks above passed despite those warnings. Those dependency
conflicts are resolved by the `.colabdeps` revision.

The fresh `.colabdeps` validation retained requests 2.32.4, Numba 0.64.0,
llvmlite 0.46.0, Bokeh 3.8.2 and NumPy 2.4.6. `pip check` reported no broken
requirements. All 962 selected installed-wheel checks passed, including LDTK's
Numba nonlinear profile calculation and serialization of EXOTIC's Bokeh image
and drawing tools. The shared binary downloads, native/Python comparisons and
joint likelihood checks again used zero C++ compiler calls and no ZIP uploads.

That setup took 73.2 seconds, compared with 91.9 seconds in the earlier run;
notebook imports took 7.4 seconds. The measured install stages were:

| Stage | Seconds |
| --- | ---: |
| Bootstrap packages (setuptools, jedi, importlib-metadata) | 21.4 |
| Dependency resolution, metadata and source preparation | 11.8 |
| Dependency wheel payload downloads | 2.9 |
| Pure Python imreg-dft wheel build | 1.4 |
| Unpacking, replacing and installing packages | 34.1 |
| Installer retrieval, shared binaries and other startup work | 1.5 |

The phase timings come from ordinary pip log boundaries. The three shared
wheels downloaded in 0.94 seconds; dependency payloads fell from approximately
187 MB to 117 MB. Package installation and bootstrap work dominated this run.
Network and VM performance vary between runtimes. Drive mounting and complete
science-image reductions were not part of this dependency validation.

Rebuild the wheel set with `scripts/build_colab_wheels.py`. It also regenerates
the exact version requirements and wheel manifest from that build's outputs.
Do not reuse an old EXOTIC wheel after changing its Python or native source;
the complete EXOTIC wheel is a snapshot of both.
