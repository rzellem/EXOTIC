# Manual validation scripts

## Claret and hybrid4 on saved light curves

`compare_hybrid4_lightcurves.py` runs both UltraNest models on an existing
`FinalLightCurve` CSV, using the stellar parameters and filter from its inits
file. The same PHOENIX profile generation supplies both intensity laws.
Run it in a foreground repository console with the existing Python environment:

```powershell
.\.venv\Scripts\python.exe scripts\build_hybrid4_cpp.py
.\.venv\Scripts\python.exe -u scripts\compare_hybrid4_lightcurves.py `
  --case 'D:\observation\inits.json' 'D:\observation\working_artifacts\FinalLightCurve_Planet_2026-04-08.csv' `
  --output 'D:\hybrid4-comparison'
```

Repeat `--case INITS_JSON FINAL_LIGHTCURVE_CSV` for multiple observations.
Only the first six numeric CSV columns are read; magnitude columns may contain
text. This is a refit of the saved, already selected and detrended photometry,
without new clipping or detrending. Both fits use the existing initial transit
bounds and duration-prior calculations. The epoch is reconstructed from the
saved BJD_TDB and phase columns, avoiding reinterpretation of legacy UTC epochs.
Recorded constant exposures and minimum live points are taken from the inits
file. If no exposure is recorded, both replay fits use instantaneous fluxes and
the limitation is stated in the report; normal reductions retain their actual
per-frame exposure arrays.

The local output directory receives paired JSON summaries, CSV data/models,
overlay PNGs and `comparison_summary.json`. Companion paths in the JSON are
relative to that report, so the complete directory can be archived or moved. The existing observations and their
outputs are not modified. Each failed case is recorded explicitly and the script
returns a non-zero exit status if any comparison failed.

## Prebuilt Colab wheels

EXOTIC wheels include the shared C library and `build_info.json`. The library
does not use the Python or NumPy C API, so its wheel is tagged `py3-none` with
the actual platform rather than a particular CPython ABI. Python 3.12 remains
the minimum supported version. An installed packaged binary is reused even
when wheel extraction gives its bundled headers newer modification times.

Build on Linux from the checkout or a genuine source distribution, using the
Python version intended for Colab. Install the build tools once in that build
environment:

```bash
python -m pip install build wheel auditwheel patchelf
python scripts/build_colab_wheels.py --output /path/to/colab-wheels
```

`--source /path/to/exotic.tar.gz` builds an existing source distribution.
The helper also collects an UltraNest wheel; use `--ultranest-version` to select
the tested version, otherwise pip selects the latest available. UltraNest's wheel is specific
to the build Python ABI. Linux compatibility tags come from `auditwheel`'s
analysis of the actual binaries. Building on an older manylinux image can
produce a wheel supporting older glibc versions; a newer Linux build must not
be relabelled as supporting an older runtime without that compatibility.
The existing `exotic-colab` notebook display helper is collected as a pure
Python wheel as well, retaining its normal package dependency.
The builder also calls `prepare_colab_dependencies.py` to collect the tested
Colab dependency wheels. It builds pure Python imreg-dft once and removes tests,
examples and HTML coverage reports from Astroquery, Photutils, colour-science
and colour-demosaicing. All retained payloads and licences are compared with
their upstream wheel files; runtime catalogues and numerical data remain.
The `1colabslim` build tag distinguishes these wheels without changing package
versions. Their regenerated RECORD protects software-installation integrity;
science data is not scanned. `dependency_audit.json` records every omission.

Transfer the resulting wheels into the Colab runtime and install their paths:

```python
from pathlib import Path
import subprocess
import sys

wheelhouse = Path("/content/colab-wheels")
wheels = [next(wheelhouse.glob("exotic-*.whl")),
          next(wheelhouse.glob("ultranest-*.whl")),
          next(wheelhouse.glob("exoticcolab-*.whl"))]
subprocess.check_call([sys.executable, "-m", "pip", "install",
                       *map(str, wheels)])
```

This uses the bundled native libraries. The source/Git installation path still
builds locally. Keep generated distributions local unless publication is
explicitly requested.

## Automatic Colab setup

`install_colab.py` installs the exact package versions advertised in a shared
wheelhouse. It uses pip's standard wheel selection and
[`--only-binary`](https://pip.pypa.io/en/stable/cli/pip_install/#cmdoption-only-binary)
for all eight bundled packages. Other dependencies continue to use
normal pip resolution. Missing downloads and incompatible native wheels are
reported directly; the installer does not retry with a Git or source build.
The installer reads the selected wheel's requirements and retains compatible
installed NumPy, requests, Numba, Bokeh, SciPy, Astropy, Matplotlib and
scikit-image versions using exact pip constraints. Missing or incompatible
versions are resolved normally. The prepared `.colabstack` wheel accepts
Colab's NumPy 2.4.6, SciPy 1.16.3, Astropy 7.2.3, Matplotlib 3.10.0 and
scikit-image 0.25.2, and preserves the original native library.
The installer passes the downloaded wheel paths directly to pip, ensuring the
prepared builds are selected. Notebook prerequisites (`jedi` and
`importlib-metadata`) join that same pip call; Step 1 has no separate bootstrap
install. Already-satisfied packages are retained without a blanket upgrade.
The prebuilt path does not need to upgrade setuptools for source builds.

The prepared wheels, version requirements and download manifest are in
[`wheelhouse/colab`](../wheelhouse/colab/README.md). The default download URLs
will work after these files are committed and pushed to the existing EXOTIC
branch. These tools do not publish or upload anything.

The complete replacement Step 1 source is in
[`colab_setup_cell.txt`](colab_setup_cell.txt). Paste it into the notebook's
existing Step 1 after the repository files have been pushed. It preserves
the Drive mount, display helpers and imports, and replaces the ZIP upload
with automatic downloads.

Once published, replace Step 1's ZIP extraction/install block with:

```python
from urllib.request import urlretrieve
from pathlib import Path
import subprocess
import sys

installer_url = (
    "https://raw.githubusercontent.com/rzellem/EXOTIC/"
    "michael_fitzgeralds_wonderful_branch_of_magic/scripts/install_colab.py"
)
installer_path = Path("/content/exotic_install_colab.py")
urlretrieve(installer_url, installer_path)
subprocess.check_call([sys.executable, str(installer_path)])
```

Keep the notebook's Drive mount, display setup and reduction cells. They do
not need a per-notebook upload or a local bundle path. Copies of the notebook
use the shared installer and current pinned wheel set. Refresh that wheel set
after changing EXOTIC's source; installing a source checkout separately would
bypass the tested binary snapshot.

For a local checkout or an already available wheelhouse:

```bash
python scripts/install_colab.py --wheelhouse wheelhouse/colab
```

`--wheelhouse` also accepts another explicitly selected shared URL. It is a
normal command argument; no application environment variables are introduced.
The build helper regenerates `requirements.txt` and `index.json` from only
the current build's repaired wheels, so previous build files are not advertised.

## C++ hybrid4 benchmark

Installation and first use build the native library automatically when needed.
`python -m exotic.api.native_build` prepares it ahead of a run; add `--force`
for an explicit rebuild. The compatibility script `build_hybrid4_cpp.py` uses
the same builder. An installed C++ compiler is required for source builds;
current builds are reused, including across concurrent bulk workers. Code,
binaries and build metadata remain in this checkout.

`benchmark_hybrid4_cpp.py` loads the shared input CSV, fixed atmosphere weights,
bounds, duration prior and exposure settings from an existing comparison report.
It runs a current Claret control and three matched-seed NumPy/C++ pairs, with
alternating backend order. Atmosphere profiles are not regenerated, and no
photometry selection, clipping, detrending or uncertainty scaling is repeated.

```powershell
.\.venv\Scripts\python.exe -u scripts\benchmark_hybrid4_cpp.py `
  --report 'D:\hybrid4-comparison\LimbDarkeningComparison_Planet_2026-04-08.json' `
  --output 'D:\hybrid4-cpp-comparison'
```

Repeat `--report` for multiple datasets. The script records per-seed posterior
intervals, goodness of fit, evidence, likelihood-call counts and full-fit timing;
it also measures flux-only and orbit/exposure-inclusive model calls. Library
loading is warmed before measuring, and compilation time is excluded. Existing
sampler settings and CPU admission behavior are retained. Numerical and
posterior differences are reported as diagnostics without new acceptance gates.
All reports, paired-model CSVs and plots stay in the requested local directory.

## C++ orbit and exposure benchmark

`benchmark_hybrid4_orbit_cpp.py` compares the previous `cpp_flux` path (Python
orbit/exposure integration, C++ flux) with `cpp` (all three calculations in C++).
Both paths use the same saved inputs and three matched seeds in alternating
order. A fresh Claret control measures the complete improvement on the same
machine; the sampler and its existing CPU admission behavior are retained.

```powershell
.\.venv\Scripts\python.exe -u scripts\benchmark_hybrid4_orbit_cpp.py `
  --report 'D:\hybrid4-comparison\LimbDarkeningComparison_Planet_2026-04-08.json' `
  --output 'D:\hybrid4-orbit-comparison'
```

Repeat `--report` for each dataset. Each local JSON contains shared atmosphere
profiles, fit settings, per-seed posteriors, reported errors, evidence, chi-square,
residual RMS, call counts and fit times. Companion CSV/PNG files compare the
first paired fits. The model comparison measures 65 identical parameter draws
and separate model-call speed; this is distinguished from complete fit speed.
Compilation is warmed before timing. Native ABI version 3 includes the Claret
quadrature model and batched orbital-posterior likelihood replay, in addition to
the hybrid4 orbit/exposure entry point, and automatically invalidates older builds.

## Full FITS reduction profiling

`benchmark_photometry_full.py` runs the same full profiling workflow with
`--backend cpp` or `--backend python`. Both controls retain the current native
transit models and orbital replay, isolating the photometry change. Use separate
output directories and run the controls sequentially so they do not compete for
CPU or disk bandwidth. Its default uses lightweight wall/stage clocks with
cProfile disabled; `--cprofile` enables per-call diagnostic profiles:

```powershell
.\.venv\Scripts\python.exe scripts\benchmark_photometry_full.py `
  --backend cpp --inits 'D:\observation\inits.json' `
  --science-dir 'D:\observation\fits' --output 'D:\comparison\cpp'
.\.venv\Scripts\python.exe scripts\benchmark_photometry_full.py `
  --backend python --inits 'D:\observation\inits.json' `
  --science-dir 'D:\observation\fits' --output 'D:\comparison\python'
```

Each output also records aperture grids/backgrounds/noise in
`aperture_grids.npz`, their star positions in `aperture_grid_index.json`, the
scout-to-sampler bounds in `bound_handoff_trace.json`, and each nested fit's
inputs and weighted posterior. Report extraction-stage timings separately from
the full elapsed time: adaptive star selection, clipping and sampler work can
change even when numerical differences are small. Full timings include the
normal network and CPU-admission waits; enabling `--cprofile` adds profiler overhead.

`profile_full_reduction.py` runs the normal complete reduction with the current
backend and a local copy of an existing initialization file. It redirects
outputs and the input-image directory, preserves the science settings, and uses
the recorded random seed or 20261121 when none was recorded. Images are read
directly; they are not copied or hashed. Supply the directory containing the
science frames directly, rather than an enclosing download/results directory.

```powershell
.\.venv\Scripts\python.exe -u scripts\profile_full_reduction.py `
  --inits 'D:\observation\inits.json' `
  --science-dir 'D:\observation\fits' `
  --output 'D:\full-reduction-profiles\Planet'
```

The foreground console remains live. Its output is also saved in `console.log`.
`profile_summary.json` records nonoverlapping main-thread stage times, function
counts, individual fits and configured sampler workers. `full_reduction.prof`
and `profile_functions.txt/.json` retain source-labeled main-thread call timings.
`profile_raw_calls.json` also keeps distinct code identities when legacy and
package imports load the same source twice; pstats source-label keys can merge
or overwrite those entries. Stage clocks remain independent of that labeling.
Callback targets absent from lsprof's top-level entries are marked as callee-only
nodes; their measured caller edges are retained and no top-level totals are invented.
Final clocks and the standard profile are saved before exporting this graph.
Regular reduction products, including the Claret/hybrid4 comparison, are written
under `reduction`. The original initialization and earlier outputs are preserved.

These are **profiled** wall times: cProfile adds overhead, and they are not a
measurement of uninstrumented production performance. Worker-thread stage calls
are counted separately; subprocess computation is not included in the parent
cProfile. High-level stage exclusive times exclude nested wrapped stages.
Inclusive function times overlap and must not be summed. External-service wait,
FITS decompression and already compiled numerical work should be identified
before deciding which Python operations would benefit from C++.

`benchmark_claret_cpp.py` separately compares the original Python and native
Claret paths on a saved comparison report's exact photometry and atmosphere
coefficients. It measures 65 identical model draws, three alternating-order
matched-seed fits, and scalar-Python/scalar-native/batched-native orbital replay
on the same positive-weight posterior draws. Replay retains the original
baseline, exposure, duration-prior and impact-parameter prior calculations and
reports score, support, weight, ESS and posterior-interval differences.

```powershell
.\.venv\Scripts\python.exe -u scripts\benchmark_claret_cpp.py `
  --report 'D:\results\LimbDarkeningComparison_Planet_date.json' `
  --output 'D:\claret-cpp-comparison'
```

Repeat `--report` for additional datasets. These saved-photometry measurements
exclude FITS reduction and cProfile overhead; compare full reductions separately.
Fit runtimes include the existing CPU admission waits; the model and replay
microbenchmarks have no admission gate. Keep system load in mind when comparing
fit wall times.

For a complete FITS reduction comparison with the same native photometry and
transit model, run `benchmark_photometry_full.py --backend cpp` once with
`--likelihood-backend python` and once with `--likelihood-backend cpp`, using
the same `--inits` and `--science-dir` and separate `--output` directories.
Both paths retain the final-plot uncertainty reuse. The default has per-call
profiling disabled; stage clocks, fit inputs/posteriors, bounds and aperture
grids are retained. Ordinary comparison JSON includes each law's actual
likelihood backend. This script runs the full reduction and does not resume
existing output automatically.

`--preparation-backend python` / `--preparation-backend cpp` compares the old
image preparation and the native histogram/repair plus lazy FITS access. Keep
`--backend cpp --likelihood-backend cpp` in both runs to isolate preparation.
Stage clocks now distinguish image reads/decompression, image headers,
calibration/demosaicing, bad-pixel repair and other frame preparation. Complete
standard outputs, aperture grids, final sampler inputs and both law summaries
remain available for scientific comparison.

## Expanded-prior warm starts

`expanded_prior_real_data_trial.py` compares an initial narrow-prior fit, an
expanded-prior warm fit, and an independent expanded-prior cold fit on saved
light curves. It retains the observation-specific priors, seeds, and comparison
criteria from the original trial. The light curves are local inputs and are
not included in this repository.

Run it from the repository root with EXOTIC installed in your environment:

```bash
python scripts/expanded_prior_real_data_trial.py \
  --dataset-path 0 "/data/NormalizedFlux_WASP-43 b_2026-03-20.txt" \
  --output /data/expanded-prior-trial.json
```

Use the index corresponding to the original observation:

| Index | Observation preset | Input format |
| --- | --- | --- |
| 0 | WASP-43b full transit, 2026-03-20 | Normalized flux |
| 1 | TrES-5b full transit, 2026-07-16 | Final light curve |
| 2 | KELT-20b partial transit, 2026-07-15 | Normalized flux |

Normalized-flux files are CSV tables with `BJD`, `Norm Flux`, `Norm Err`, and
`AM` headers. Final-light-curve files are CSV tables with time in column 1,
flux in column 3, flux error in column 4, and airmass in column 6; comment lines
start with `#`.

Repeat `--dataset-path INDEX PATH` to supply multiple datasets. By default,
the script runs the supplied indices in argument order. The existing
`--dataset-index INDEX` option can select a subset of those inputs and can
also be repeated. A path is required for every selected dataset. Relative
paths are relative to the directory from which you run the command.

`--live-points` and `--seed` retain their existing defaults of 100 and 4729.
The JSON output retains the fit summaries, comparison diagnostics, and
pass/fail result. Loading or fitting failures are reported directly.
