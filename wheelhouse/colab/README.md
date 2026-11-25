# Prebuilt Colab software

The wheelhouse pins a build of EXOTIC and the dependency set validated in the
Colab CPU runtime on 22 November 2026. The validation results below describe
that November 22 setup. The current EXOTIC source revision is in the wheel's
version, recorded with the dependency versions and filenames in
`requirements.txt` and `index.json`.

The complete set supports CPython 3.13, Linux x86_64 and glibc 2.38 or newer.
EXOTIC's ctypes library uses native ABI 6 and has no Python C ABI dependency.
UltraNest's extension uses CPython 3.13. The original wheel license files are
preserved.

The current manifest includes 22 wheels: the three EXOTIC/UltraNest/display
packages, pure Python imreg-dft 2.0.0, and trimmed upstream Astroquery 0.4.11,
Photutils 3.1.1, colour-science 0.4.7 and colour-demosaicing 0.2.7 wheels. Their
`1colabslim` build tag preserves the upstream versions. Only tests, examples
and HTML coverage reports are omitted; every retained payload and licence is
unchanged, except for the WHEEL build tag and regenerated RECORD. Catalogues,
runtime numerical data and compiled extensions remain. RECORD is regenerated
once for software-distribution integrity, with no science-file hashing.
Version/build identifiers cannot replace the wheel format's required member
list and per-member checksums after those development files are omitted.
`dependency_audit.json` records the omitted files and their sizes.

The other 14 wheels are unmodified upstream Jedi, Astroalign, barycorrpy,
Bottleneck, LDTk, lmfit, asteval, psutil, PyVO, REBOUND, SEP, uncertainties,
corner and jplephem, pinned to the versions used in the preceding successful
fresh-runtime validation. Shared Colab libraries are supplied by the runtime.
Four workers download the manifest wheels concurrently, with stable report
ordering. The installer records download wall time separately from overlapping
per-wheel times. It then uses one pip call with `--no-index --no-compile`:
dependency compatibility checks remain enabled, while PyPI lookups and
installation-time Python bytecode compilation are avoided. Python may prepare
bytecode during first imports, so compare installation plus imports.

The fresh 22-wheel validation downloaded all 34.0 MB in 1.31 seconds of wall
time, using four workers, and made no PyPI package downloads. Pip resolved
dependencies locally, kept all eight protected Colab versions, and `pip check`
reported no broken requirements. No bytecode files existed in the checked
newly installed packages before their ordinary imports. All 1,209 selected
installed-wheel checks passed, with zero source builds, C++ compiler calls or
ZIP uploads. Local installer/dependency/native checks passed 32 tests, with
two Colab-only checks skipped locally and exercised in that fresh runtime.

| Software startup stage | Previous 8-wheel setup | Complete 22-wheel setup |
| --- | ---: | ---: |
| Installation | 20.50 s | 14.67 s |
| First imports | 7.88 s | 10.64 s |
| Installation plus first imports | 28.38 s | 25.30 s |

Installation was 28.4% shorter, while the combined install/import measurement
was 10.8% shorter (3.08 seconds saved). First imports were longer, so the
installation-only improvement overstates the observed complete startup gain.
The three changes were measured together; their individual effects were not
isolated. These are separate fresh-VM observations, excluding VM allocation,
Drive mounting and validation-only reference installs/refits. Pip itself took
12.65 seconds: approximately 6.47 seconds for startup/local dependency
resolution and 6.18 seconds for unpack/install.

The new setup completed 24 matched HAT-P-26b/WASP-16b Claret/hybrid4 refits,
using three seeds and the same stored photometry, atmosphere coefficients,
exposures, masks, priors and sampling settings. All 12 pairs against isolated
untrimmed dependencies had identical complete posterior summaries, reported
errors, maximum-likelihood/reported parameters, models, RMS, chi-square,
evidence and its error, bounds and call counts. The 1,028 likelihood/prior
replays and Bayer/imreg dependency probes also matched exactly. No fit reached
its call limit. The same fields, models and prior bounds also matched all 12
saved fits from the preceding eight-wheel setup exactly. Reference installs
and refits are excluded from startup timing; images and atmosphere profiles
were not regenerated.

| Package | Original wheel MB | Prepared wheel MB |
| --- | ---: | ---: |
| Astroquery | 11.139 | 0.741 |
| Photutils | 10.075 | 9.277 |
| colour-science | 9.066 | 1.644 |
| colour-demosaicing | 3.263 | 0.015 |

Together these omit 21.866 MB compressed and 107.803 MB unpacked, across 1,786
development files. The pure Python imreg-dft wheel is 47 KB and avoids its
per-runtime source preparation/build. The installer supplies these exact files
and the notebook's Jedi/importlib-metadata requirements to one pip call,
retaining compatible installed packages without a blanket upgrade. The normal
Step 1 no longer has a separate bootstrap pip command. Fresh-runtime timing
and installed-wheel validation used the pushed eight-wheel manifest directly.

Normal installation took 20.5 seconds versus 42.4 seconds in the preceding
fresh VM: 21.9 seconds, or 51.6%, shorter in this pair. The single pip call took
17.4 seconds; log-boundary phases were 7.5 seconds for startup/resolution,
3.4 seconds for PyPI payload downloads, zero for source builds, and 6.5 seconds
for unpack/install. The eight shared wheels downloaded in 2.3 seconds, with
the remaining time in manifest/installer retrieval and startup. Imports took
7.9 seconds separately. Allocation, Drive mounting and validation-only work
are excluded; separate-VM and network performance vary.

Measured wheel downloads totalled 34.0 MB: 15.1 MB from the shared wheelhouse
and 18.9 MB from PyPI, including the 4.9 MB Jedi wheel. The previous quiet
bootstrap did not record its payload bytes, so its 47.6 MB dependency figure
is not directly comparable to this combined total. The four trimmed wheels'
21.866 MB reduction above compares the actual upstream/prepared files.

All 1,209 selected installed-wheel checks passed with zero source builds,
C++ compiler calls or ZIP uploads. These include all three Bayer demosaicing
algorithms and imreg-dft translation, alongside the native/Python models,
photometry, plotting, joint fits, timing and FITS/WCS checks. `pip check`
reported no broken requirements, and all eight protected Colab versions were
retained. The local full suite passed 3,108 tests with two Colab/IPython checks
skipped; four additional image API checks passed separately.

The scientific comparison completed 24 independent fits of recorded HAT-P-26b
(241 points) and WASP-16b (72 points) photometry: Claret and hybrid4 under both
slim and isolated untrimmed copies of the same dependency versions, using
three matched seeds. Atmosphere coefficients, photometry, masks, exposures,
priors and normal sampling controls were fixed. All 12 pairs had identical
complete posterior summaries, maximum-likelihood/reported parameters,
reported errors and adjustments, models, RMS, chi-square, evidence and its
error, prior bounds, likelihood-call counts and sampling budgets. No fit
reached its call limit. The 1,028 prior/likelihood replay points, three
demosaicing probes and image-translation probe also matched exactly. The normal runtime
kept its prepared libraries and `pip check` remained clean. This comparison
used saved photometry and did not regenerate profiles or reduce images.
Reference downloads and all refits were excluded from installation timing.

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
requirements. Newer compatible versions remain supported.

The fresh `.colabstack` runtime retained all eight requested packages and their
shared dependencies; `pip check` reported no broken requirements. All 1,205
selected installed-wheel checks passed, including plotting, image alignment
and downsampling, FITS/WCS and timing, LDTK, joint fits and native/reference
maths. The local suite with Astropy 7.2.3, Matplotlib 3.10.0, scikit-image
0.25.2 and SciPy 1.16.3 passed 3,097 tests; two Colab/IPython checks were skipped
locally and passed in Colab.

Normal setup took 42.4 seconds versus 55.7 seconds in the previous fresh
runtime (13.3 seconds, or 24%, shorter in this pair). Dependency payloads fell
from 81.6 MB to 47.6 MB, avoiding 34.0 MB of Astropy/Matplotlib/scikit-image
downloads. Pip installation fell from 37.9 to 25.8 seconds and its
unpack/replace/install phase from 23.2 to 9.4 seconds. Bootstrap packages took
15.1 seconds; notebook imports took 9.0 seconds. The three shared binary wheels
downloaded in 0.64 seconds, without native compilation or ZIP uploads. These
are separate fresh-VM measurements; runtime and network performance vary.
VM allocation, Drive mounting, validation-only reference installs and
scientific refits are excluded from normal setup timing.

Scientific validation performed 24 independent refits: HAT-P-26b (241 points)
and WASP-16b (72 points), each with Claret and hybrid4, using three matched
seeds under both the retained stack and an isolated reference stack of
Astropy 8.0.3, Matplotlib 3.11.2 and scikit-image 0.26.0. NumPy 2.4.6 and SciPy
1.16.3 were shared. All 12 pairs had identical complete posterior summaries,
maximum-likelihood and reported parameters, reported errors and adjustments,
models, RMS, chi-square, evidence and evidence errors, prior bounds,
likelihood-call counts and sampling budgets. No fit reached its call limit.
The 1,028 prior/likelihood replay points also matched exactly. Photometry,
stored atmosphere coefficients, masks, exposure integration and priors were
held fixed; no images were re-reduced or atmosphere profiles regenerated.
The normal runtime retained its original libraries and `pip check` remained
clean after validation. The notebook's diagnostic was restored to its small
optional model check; normal Run all does not repeat these benchmark fits.

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
