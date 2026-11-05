"""
EXOTIC GPU — Installation & Environment Validator (Cross-platform)
===================================================================
Usage (any OS: Linux / Windows / macOS):

    # From the PARENT folder of EXOTIC/ (official repository root):
    python EXOTIC/check_gpu_install.py

    # Or from inside EXOTIC/:
    python check_gpu_install.py

Note: this script does NOT depend on any concrete distribution, hard-coded
path, or specific GPU model. It generically probes the environment and
offers installation alternatives for every platform.

Validates step by step:
  1. Python version (3.10–3.12 recommended, 3.9–3.14 best-effort)
  2. Core packages (numpy, scipy, astropy, photutils, scikit-image, …)
  3. EXOTIC imports correctly
  4. Available GPU backends (PyTorch CUDA / CuPy / Numba CUDA)
  5. EXOTIC.exotic_gpu loads and reports status
  6. Quick phase_cross_correlation benchmark on 512×512
  7. Quick compute_star_aperture_grid benchmark

Exit code semantics: 0 = environment OK, 1 = missing steps.
"""

import sys
import os
import time
from types import SimpleNamespace

RESULTS = []


# -----------------------------------------------------------------
# Cross-platform helpers (no hardcodes, no OS-conditional branching)
# -----------------------------------------------------------------
def section(title):
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def ok(msg):
    print(f"  OK  {msg}")
    RESULTS.append(("OK", msg))


def warn(msg):
    print(f"  !!  {msg}")
    RESULTS.append(("WARN", msg))


def fail(msg):
    print(f"  FAIL {msg}")
    RESULTS.append(("FAIL", msg))


def info(msg):
    print(f"  ·   {msg}")


def env_var_set(name: str) -> str:
    """Return the shell line that sets an env var according to the OS / shell."""
    if os.name == "nt":
        return f"set {name}=1"
    sh = os.environ.get("SHELL", "").lower()
    if "fish" in sh:
        return f"set -gx {name} 1"
    return f"export {name}=1"


def env_var_unset(name: str) -> str:
    if os.name == "nt":
        return f"set {name}=0"
    sh = os.environ.get("SHELL", "").lower()
    if "fish" in sh:
        return f"set -gx {name} 0"
    return f"export {name}=0"


# -----------------------------------------------------------------
# 1. Python version
# -----------------------------------------------------------------
section("1. Python version")
# Supported range: 3.9 – 3.14 (strictly recommended: 3.10 – 3.12 per EXOTIC)
PY_STRICT = (3, 10) <= sys.version_info <= (3, 12)
PY_LOOSE = (3, 9) <= sys.version_info < (3, 15)
info(f"Platform       : {sys.platform}  ({os.name})")
info(f"Python version : {sys.version.split()[0]}")
info(f"Executable     : {sys.executable}")
if PY_STRICT:
    ok(f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro} — recommended range")
elif PY_LOOSE:
    warn(f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro} works but is outside the recommended range (3.10–3.12).")
else:
    if sys.version_info < (3, 9):
        fail("Python too old (< 3.9). EXOTIC and some deps will not build. Install 3.10.x or 3.11.x.")
    else:
        warn(f"Python 3.{sys.version_info.minor} is very new. Some pinned deps (e.g. numpy 1.26.x) may not be available yet.")

# -----------------------------------------------------------------
# 2. Core required packages
# -----------------------------------------------------------------
section("2. EXOTIC core dependencies")
REQ_CORE = [
    ("numpy", ">=1.24 recommended"),
    ("scipy", ">=1.11"),
    ("astropy", ">=6.0"),
    ("photutils", ">=1.10"),
    ("skimage", ">=0.22 (package scikit-image)"),
    ("matplotlib", ">=3.8"),
    ("numba", "JIT CPU acceleration"),
]
for pkg, _ver in REQ_CORE:
    try:
        mod = __import__(pkg)
        ver = getattr(mod, "__version__", "?")
        ok(f"{pkg:15s} {ver}")
    except ImportError:
        fail(f"{pkg:15s} NOT INSTALLED → pip install {pkg}")
        if pkg == "skimage":
            info("    (real PyPI name is scikit-image, imported as skimage)")

# -----------------------------------------------------------------
# 3. Available GPU backends (cross-platform)
# -----------------------------------------------------------------
section("3. Available GPU backends")
GPU_BACKENDS = {}
try:
    import torch
    GPU_BACKENDS["pytorch"] = {
        "installed": True,
        "version": getattr(torch, "__version__", "?"),
        "cuda": bool(getattr(torch.cuda, "is_available", lambda: False)()),
        "device_name": None,
        "memory_gb": 0.0,
    }
    if GPU_BACKENDS["pytorch"]["cuda"]:
        try:
            GPU_BACKENDS["pytorch"]["device_name"] = torch.cuda.get_device_name(0)
            GPU_BACKENDS["pytorch"]["memory_gb"] = round(
                torch.cuda.get_device_properties(0).total_memory / 1e9, 2
            )
        except Exception:
            pass
        ok(
            f"PyTorch {GPU_BACKENDS['pytorch']['version']} — CUDA: "
            f"{GPU_BACKENDS['pytorch']['device_name'] or 'GPU'} "
            f"({GPU_BACKENDS['pytorch']['memory_gb']} GB VRAM)"
        )
    else:
        warn(
            f"PyTorch {GPU_BACKENDS['pytorch']['version']} installed but CUDA NOT detected."
        )
        info("  Install with CUDA support (any Linux/Windows):")
        info("    pip install torch --index-url https://download.pytorch.org/whl/cu121")
        info("  (replace cu121 by cu118 if you are pinned to legacy CUDA 11.8 drivers)")
except ImportError:
    GPU_BACKENDS["pytorch"] = {"installed": False}
    warn("PyTorch CUDA NOT installed. (RECOMMENDED for FFT alignment and batch photometry.)")
    info("    pip install torch --index-url https://download.pytorch.org/whl/cu121")

try:
    import cupy
    mem_info = None
    device_name = None
    try:
        mem_info = cupy.cuda.Device().mem_info
        device_name = cupy.cuda.Device().name
    except Exception:
        mem_info = None
        device_name = None
    GPU_BACKENDS["cupy"] = {
        "installed": True,
        "cuda": mem_info is not None,
        "version": getattr(cupy, "__version__", "?"),
        "device_name": device_name,
        "memory_gb": round(mem_info[1] / 1e9, 2) if mem_info else 0.0,
    }
    if GPU_BACKENDS["cupy"]["cuda"]:
        ok(
            f"CuPy {GPU_BACKENDS['cupy']['version']} — CUDA: "
            f"{GPU_BACKENDS['cupy']['device_name'] or 'GPU'} "
            f"({GPU_BACKENDS['cupy']['memory_gb']} GB VRAM)"
        )
    else:
        warn(f"CuPy {GPU_BACKENDS['cupy']['version']} installed but no CUDA GPU detected.")
except ImportError:
    GPU_BACKENDS["cupy"] = {"installed": False}
    info("CuPy CUDA (alternative to PyTorch for even faster FFTs):")
    info("    pip install cupy-cuda12x   # (CUDA 12.x)")
    info("    pip install cupy-cuda11x   # (CUDA 11.x)")

try:
    from numba import cuda as _nb_cuda
    GPU_BACKENDS["numba"] = {
        "installed": True,
        "cuda": bool(_nb_cuda.is_available()),
        "version": None,
    }
    try:
        import numba
        GPU_BACKENDS["numba"]["version"] = numba.__version__
    except Exception:
        pass
    if GPU_BACKENDS["numba"]["cuda"]:
        ok(f"Numba CUDA {GPU_BACKENDS['numba']['version'] or ''} detected")
    else:
        info("Numba available (bundled JIT CPU acceleration already used by EXOTIC base)")
except Exception:
    GPU_BACKENDS["numba"] = {"installed": False}

HAS_ANY_GPU = any(v.get("cuda") for v in GPU_BACKENDS.values() if v.get("installed"))
if not HAS_ANY_GPU:
    info("No CUDA backend detected. EXOTIC will run on CPU;")
    info("GPU code is automatically disabled (robust CPU fallback).")

# -----------------------------------------------------------------
# 4. Import EXOTIC
# -----------------------------------------------------------------
section("4. Import EXOTIC main module")

# Ensure both EXOTIC/ root and its parent are on sys.path for every possible
# way the user can invoke this script.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR = os.path.dirname(_THIS_DIR)
for _p in (_ROOT_DIR, _THIS_DIR):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)

EXOTIC_IMPORTED = False
EXOTIC_MODULE = None
try:
    import EXOTIC.exotic.exotic as exo_mod
    EXOTIC_IMPORTED = True
    EXOTIC_MODULE = exo_mod
    ok("EXOTIC.exotic.exotic imported successfully")
    info(f"  File: {exo_mod.__file__}")
    # Check GPU injection variables
    if hasattr(exo_mod, "_EXOTIC_GPU_AVAILABLE"):
        ok("_EXOTIC_GPU_AVAILABLE marker present in exotic.py")
        info(f"  _EXOTIC_GPU_AVAILABLE    = {exo_mod._EXOTIC_GPU_AVAILABLE}")
        info(f"  _EXOTIC_GPU_ALIGN_BACKEND = {exo_mod._EXOTIC_GPU_ALIGN_BACKEND}")
        info(f"  _EXOTIC_GPU_PHOT_BACKEND  = {exo_mod._EXOTIC_GPU_PHOT_BACKEND}")
    else:
        warn("GPU injection block missing from exotic.py. Make sure you are on the GPU-acceleration branch.")
except Exception as e:
    try:
        # Alternate layout: top-level exotic package (no EXOTIC. prefix)
        import exotic.exotic as exo_mod2
        EXOTIC_IMPORTED = True
        EXOTIC_MODULE = exo_mod2
        ok("exotic.exotic imported (top-level installed package, no EXOTIC. prefix)")
    except Exception as e2:
        fail(f"Failed to import EXOTIC: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc(limit=3)
        info("  Make sure you run this script FROM the PARENT folder of EXOTIC/:")
        info("    cd /path/to/EXOTIC-GPU   # (the folder that contains EXOTIC/)")
        info("    python EXOTIC/check_gpu_install.py")

# -----------------------------------------------------------------
# 5. Import EXOTIC.exotic_gpu
# -----------------------------------------------------------------
section("5. Import EXOTIC.exotic_gpu and probe enable_exotic_gpu()")
EXOTIC_GPU_IMPORTED = False
STATUS = None
if not EXOTIC_IMPORTED:
    warn("Skipping (EXOTIC main module did not import cleanly)")
else:
    try:
        from EXOTIC.exotic_gpu import enable_exotic_gpu, get_full_gpu_status
        EXOTIC_GPU_IMPORTED = True
        st = get_full_gpu_status()
        ok("EXOTIC.exotic_gpu imported successfully")
        info(f"  alignment backend : {st['alignment']['backend']} "
             f"(available={st['alignment']['available']})")
        info(f"  photometry backend: {st['photometry_backend']} (gpu={st['photometry_gpu']})")
        info(f"  psf backend       : {st['psf_backend']} (gpu={st['psf_gpu']})")
        # Activate
        STATUS = enable_exotic_gpu(verbose=False, force=(not HAS_ANY_GPU))
        if STATUS.get("enabled"):
            ok("enable_exotic_gpu() returned enabled=True")
        else:
            if HAS_ANY_GPU:
                warn("enable_exotic_gpu() returned enabled=False even though a GPU was detected. Check installation.")
            else:
                info("enable_exotic_gpu() returned enabled=False (no GPU detected — CPU-optimised mode only).")
    except Exception as e:
        fail(f"Failed to import EXOTIC.exotic_gpu: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc(limit=4)

# -----------------------------------------------------------------
# 6. Benchmark: Alignment phase_cross_correlation (1 img 512x512, N=10 iters)
# -----------------------------------------------------------------
section("6. Quick benchmark: phase_cross_correlation alignment")
BENCH_SHAPE = (512, 512)
N_ALIGN = 10
if not EXOTIC_IMPORTED or not hasattr(EXOTIC_MODULE, "phase_cross_correlation"):
    warn("Skipping (EXOTIC.phase_cross_correlation unavailable)")
else:
    import numpy as np
    try:
        rng = np.random.default_rng(42)
        ref = rng.normal(1000, 20, BENCH_SHAPE).astype(np.float64)
        mov = np.roll(np.roll(ref, 17, axis=0), -23, axis=1)

        # Warmup
        EXOTIC_MODULE.phase_cross_correlation(mov, ref, 1)
        t0 = time.perf_counter()
        for _ in range(N_ALIGN):
            sh, err, _ = EXOTIC_MODULE.phase_cross_correlation(mov, ref, 1)
        dt = time.perf_counter() - t0
        sh = np.asarray(sh).flatten()
        ok(f"phase_cross_correlation (512×512): {dt:.3f} s total / {dt/N_ALIGN*1000:.1f} ms/call")
        info(f"  Measured shift: {sh[:2].tolist()}  Expected: [±17, ∓23]  Error={err:.4f}")
        drift = abs(abs(sh[0]) - 17) + abs(abs(sh[1]) - 23)
        if drift < 0.5:
            ok("Numerical result correct (drift < 0.5 px)")
        else:
            warn(f"Shift drift = {drift:.2f} px (> 0.5). Check upsample_factor or FFT backend.")
    except Exception as e:
        fail(f"Alignment benchmark failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc(limit=4)

# -----------------------------------------------------------------
# 7. Benchmark: Aperture photometry (1 star, 6ap×3an × 200 calls)
# -----------------------------------------------------------------
section("7. Quick benchmark: compute_star_aperture_grid (1 star)")
if not EXOTIC_IMPORTED or not hasattr(EXOTIC_MODULE, "compute_star_aperture_grid"):
    warn("Skipping (EXOTIC.compute_star_aperture_grid unavailable)")
else:
    import numpy as np
    try:
        rng = np.random.default_rng(7)
        img = rng.normal(1000, 20, (400, 400)).astype(np.float32)
        # Synthetic star
        yy, xx = np.ogrid[:400, :400]
        xc, yc = 213.7, 189.3
        r2 = (xx - xc) ** 2 + (yy - yc) ** 2
        img += 20000 * np.exp(-r2 / (2 * 2.1 ** 2))
        aps = np.array([1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
        ans = np.array([3.0, 5.0, 8.0])

        # Warmup
        EXOTIC_MODULE.compute_star_aperture_grid(img, 0, xc, yc, aps, ans, fast_mode=True)
        t0 = time.perf_counter()
        N_CALLS = 200
        for _ in range(N_CALLS):
            fg, bg = EXOTIC_MODULE.compute_star_aperture_grid(img, 0, xc, yc, aps, ans, fast_mode=True)
        dt = time.perf_counter() - t0
        ok(f"compute_star_aperture_grid (6ap×3an, fast_mode): {dt:.3f} s total / {dt/N_CALLS*1000:.2f} ms/call")
        info(f"  flux_grid[3,1]={fg[3,1]:.1f}  bg_grid[3,1]={bg[3,1]:.1f}")
        if np.all(np.isfinite(fg)):
            ok("Flux grid numerically valid (all finite)")
        else:
            warn(f"Found {np.sum(~np.isfinite(fg))} non-finite entries in flux_grid.")
    except Exception as e:
        fail(f"Photometry benchmark failed: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc(limit=4)

# -----------------------------------------------------------------
# 8. Summary
# -----------------------------------------------------------------
section("8. ENVIRONMENT SUMMARY")
N_OK = sum(1 for x, _ in RESULTS if x == "OK")
N_WARN = sum(1 for x, _ in RESULTS if x == "WARN")
N_FAIL = sum(1 for x, _ in RESULTS if x == "FAIL")
print(f"  OKs  : {N_OK}")
print(f"  WARN : {N_WARN}")
print(f"  FAIL : {N_FAIL}")
print()
if N_FAIL == 0:
    print("  ENVIRONMENT READY TO RUN EXOTIC")
    print()
    print("  To run the full pipeline WITH GPU acceleration:")
    print(f"    {env_var_set('EXOTIC_ENABLE_GPU')}")
    print(f"    python -m EXOTIC.exotic.exotic --inits your_inits.json")
    print()
    print("  To disable GPU and get the CPU baseline (speedup comparison):")
    print(f"    {env_var_unset('EXOTIC_ENABLE_GPU')}")
    print(f"    python -m EXOTIC.exotic.exotic --inits your_inits.json")
    print()
    print("  (Optional) Manually toggle the GPU patch from your script:")
    print("    from EXOTIC.exotic_gpu import enable_exotic_gpu")
    print("    enable_exotic_gpu()   # prints the GPU backend banner")
else:
    print("  FAILURES PRESENT — fix the FAIL lines above before running EXOTIC.")
    print()
    print("  Recommended dependency installation path (all OSes):")
    print("    # 1) EXOTIC base deps:")
    print("    pip install -r EXOTIC/requirements.txt")
    print()
    print("    # 2) Recommended GPU backend (pick ONE):")
    print("    pip install torch --index-url https://download.pytorch.org/whl/cu121")
    print("    # Or:")
    print("    pip install cupy-cuda12x")
    print()
    if os.name != "nt":
        print("  Linux hint: verify the NVIDIA driver first with:")
        print("    nvidia-smi   # must show your GPU card name + VRAM + driver version")
        print("    • Stock Arch / Manjaro / EndeavourOS (DKMS kernels):")
        print("        sudo pacman -S nvidia-dkms cuda")
        print()
        print("    • CachyOS / kernels with BAKED-IN NVIDIA modules:")
        print("        !! DO NOT INSTALL nvidia-dkms — it conflicts with")
        print("           linux-cachyos-*-nvidia-open / -nvidia packages that")
        print("           already ship the precompiled module.")
        print("        If you need the CUDA toolkit only:  sudo pacman -S cuda")
        print()
        print("    • Ubuntu 24.04 / Debian 12 / Mint:")
        print("        sudo apt install nvidia-driver-560 nvidia-cuda-toolkit")

sys.exit(0 if N_FAIL == 0 else 1)
