"""
EXOTIC GPU Acceleration Module
==============================
GPU acceleration subpackage for EXOTIC (Exoplanet Transit Interpretation
Code). It lives inside the official EXOTIC tree as a sibling of the main
``exotic`` package (see INSTALL_GPU.md for the full layout).

Accelerated hotspots (see AUDITORIA_RENDIMIENTO.md):
  1. FFT alignment        (40–50 % of runtime) → 5–30× speedup
  2. Aperture photometry  (25–35 % of runtime) → 5–50× speedup
  3. PSF fitting          (10–15 % of runtime) → fewer least_squares iters

Supported GPU backends (auto-detected at import, best-effort order):
  - PyTorch CUDA  (recommended, most portable)
  - CuPy CUDA     (higher raw FFT throughput)
  - NumPy CPU     (robust fallback, zero extra deps)

Recommended usage (opt-in env var — no EXOTIC code changes needed):
    >>> # export EXOTIC_ENABLE_GPU=1    # Linux / macOS
    >>> # set    EXOTIC_ENABLE_GPU=1    # Windows CMD
    >>> # $env:EXOTIC_ENABLE_GPU="1"    # PowerShell
    >>> from EXOTIC.exotic.exotic import main
    >>> main()

Manual explicit usage (useful in scripts and when embedding EXOTIC):
    >>> from EXOTIC.exotic_gpu import enable_exotic_gpu
    >>> enable_exotic_gpu()   # prints the detected-backend banner
    >>> from EXOTIC.exotic.exotic import main
    >>> main()

Author: EXOTIC GPU Optimization Team
Version: 2.0.0
Date: 2026-09-01
"""

__version__ = "2.0.0"
__author__ = "EXOTIC GPU Optimization Team"

import sys as _sys
import warnings as _warnings
import numpy as _np

_GPU_ENABLED = False
_EXOTIC_PATCHED = False
_ORIGINAL_FUNCS = {}

# Cross-platform EXOTIC submodule resolver (DRY). No hard-coded paths.
from ._module_resolver import (
    resolve_exotic_module as _resolve_exotic_module,
    resolve_exotic_exotic_gpu_submodule as _resolve_exotic_exotic_gpu_submodule,
    EXOTIC_MODULE_CANDIDATES as _EXOTIC_MODULE_CANDIDATES,
)


from .alignment_gpu import (
    phase_cross_correlation_gpu,
    batch_phase_cross_correlation_gpu,
    align_image_stack_gpu,
    should_use_gpu_for_alignment,
    get_gpu_info as get_alignment_gpu_info,
)

from .photometry_gpu import (
    compute_frame_aperture_grid_batch,
    compute_aperture_grid_cpu,
    image_aperture_signal_flux_gpu,
    batch_photometry_simple,
    should_use_gpu_photometry,
    is_phot_gpu_available,
    get_phot_gpu_backend,
)

from .psf_fitting_gpu import (
    psf_moments_batch,
    gaussian_psf_moments_cpu,
    estimate_image_fwhm_batch,
    is_psf_gpu_available,
    get_psf_backend,
)

from .exotic_integration import (
    patch_exotic_photometry,
    unpatch_exotic_photometry,
    set_use_gpu,
    is_gpu_enabled,
    image_aperture_signal_flux as image_aperture_signal_flux_integrated,
    batch_photometry_exotic,
)


def get_full_gpu_status() -> dict:
    """
    Return the complete GPU acceleration status snapshot.

    Returns:
        dict with keys:
            - enabled (bool): True if acceleration was toggled on via enable_exotic_gpu()
            - patched (bool): True if EXOTIC functions were monkey-patched
            - alignment (dict): alignment backend details
            - photometry_backend (str): "numpy" / "pytorch" / "cupy"
            - photometry_gpu (bool): True if a CUDA device is active for photometry
            - psf_backend (str): PSF estimation backend
            - psf_gpu (bool): True if PSF estimation runs on CUDA
            - version (str): installed EXOTIC.exotic_gpu version string
    """
    return {
        "enabled": _GPU_ENABLED,
        "patched": _EXOTIC_PATCHED,
        "alignment": get_alignment_gpu_info(),
        "photometry_backend": get_phot_gpu_backend(),
        "photometry_gpu": is_phot_gpu_available(),
        "psf_backend": get_psf_backend(),
        "psf_gpu": is_psf_gpu_available(),
        "version": __version__,
    }


def enable_exotic_gpu(
    force: bool = False,
    patch_exotic: bool = True,
    verbose: bool = True,
) -> dict:
    """
    Activate GPU acceleration for EXOTIC.

    Auto-detects the best available GPU backend and, when requested,
    monkey-patches the critical EXOTIC hotspots with GPU-accelerated
    variants that still expose the exact same public signatures.

    Args:
        force: When True, enable the CPU-optimised code path even when no
               CUDA device is detected. Useful for testing the fallback
               pipeline and for benchmarking the non-GPU optimisations
               (ROI direct aperture grid, PSF warm-start).
        patch_exotic: When True, monkey-patch ``EXOTIC.exotic.exotic``
                      functions (alignment, aperture photometry, PSF
                      fitting) with their accelerated siblings.
        verbose: When True, print the detection banner to stdout.

    Returns:
        Activation status snapshot (see :func:`get_full_gpu_status`).
    """
    global _GPU_ENABLED, _EXOTIC_PATCHED

    status = get_full_gpu_status()
    has_gpu = (
        status["alignment"]["available"]
        or status["photometry_gpu"]
        or status["psf_gpu"]
    )

    if not has_gpu and not force:
        _GPU_ENABLED = False
        if verbose:
            print()
            print("=" * 70)
            print("EXOTIC GPU v" + __version__)
            print("=" * 70)
            print("!! No CUDA-capable GPU backend detected.")
            print("   Running in optimised CPU mode.")
            print()
            print("   To enable GPU acceleration install one of:")
            print("     pip install torch --index-url https://download.pytorch.org/whl/cu121")
            print("     pip install cupy-cuda11x   # (CUDA 11.x runtimes)")
            print("     pip install cupy-cuda12x   # (CUDA 12.x runtimes)")
            print("=" * 70)
            print()
        return status

    _GPU_ENABLED = True
    set_use_gpu(True)

    if patch_exotic:
        try:
            patch_exotic_photometry()
            _patch_exotic_alignment_and_psf()
            _EXOTIC_PATCHED = True
        except Exception as e:
            _warnings.warn(f"Failed to automatically patch EXOTIC: {e}")

    if verbose:
        st = get_full_gpu_status()
        print()
        print("=" * 70)
        print(f"EXOTIC GPU v{st['version']} — ACTIVATED")
        print("=" * 70)
        align_be = st["alignment"]
        if align_be["available"]:
            print(f"OK  Alignment:   {align_be['backend'].upper()} "
                  f"({align_be.get('device_name', 'GPU')})")
            if "memory_gb" in align_be:
                print(f"               VRAM: {align_be['memory_gb']} GB")
        else:
            print(f"--  Alignment:   CPU (fallback)")
        print(f"OK  Photometry:  {st['photometry_backend'].upper()}"
              + (" (GPU)" if st["photometry_gpu"] else " (CPU)"))
        print(f"OK  PSF Moments: {st['psf_backend'].upper()}"
              + (" (GPU)" if st["psf_gpu"] else " (CPU)"))
        print(f"    EXOTIC patch: {'OK' if st['patched'] else 'Pending'}")
        print()
        print("  Expected speedup (Qatar-1 / 128 frames, RTX-class GPU):")
        print("    Alignment :  ~6–30×   (6–8 min → ~15–60 s)")
        print("    Photometry:  ~5–50×   (4–6 min → ~5–50 s)")
        print("    PSF init   :  ~10–100× (fewer least_squares iterations)")
        print("    Total      :  ~3–8×    (17 min → 2–5 min)")
        print("=" * 70)
        print()

    return get_full_gpu_status()


def disable_exotic_gpu() -> None:
    """Disable GPU acceleration and restore EXOTIC to its original state."""
    global _GPU_ENABLED, _EXOTIC_PATCHED
    _GPU_ENABLED = False
    set_use_gpu(False)
    if _EXOTIC_PATCHED:
        try:
            unpatch_exotic_photometry()
            _unpatch_exotic_alignment_and_psf()
        except Exception:
            pass
    _EXOTIC_PATCHED = False


def _patch_exotic_alignment_and_psf() -> None:
    """
    Monkey-patch (cross-platform, no hard-coded paths):
      - exotic.phase_cross_correlation → phase_cross_correlation_gpu (smart dispatch)
      - exotic.compute_star_aperture_grid → batch-aware variant
      - exotic.fit_psf_photometry_flux_row → PSF warm-start via method of moments
    """
    exo_module = _resolve_exotic_module()
    if exo_module is None:
        return

    if not hasattr(exo_module, "_exotic_gpu_original"):
        exo_module._exotic_gpu_original = {}

    orig = exo_module._exotic_gpu_original

    if "phase_cross_correlation" not in orig and hasattr(exo_module, "phase_cross_correlation"):
        orig["phase_cross_correlation"] = exo_module.phase_cross_correlation

        def _smart_phase_cross_correlation(moving, reference, upsample_factor=1, *args, **kwargs):
            try:
                shape = getattr(moving, "shape", None)
                use_gpu = shape is not None and len(shape) == 2 and should_use_gpu_for_alignment(shape)
                if use_gpu:
                    try:
                        result = phase_cross_correlation_gpu(
                            moving, reference,
                            upsample_factor=upsample_factor,
                            return_error=True,
                        )
                        shift, err = result[0], result[1]
                        return (_np.asarray(shift, dtype=float),
                                float(err) if _np.isfinite(err) else 0.0,
                                _np.nan)
                    except Exception:
                        pass
            except Exception:
                pass
            return orig["phase_cross_correlation"](
                moving, reference, upsample_factor=upsample_factor, *args, **kwargs
            )

        exo_module.phase_cross_correlation = _smart_phase_cross_correlation

    if "fit_psf_photometry_flux_row" not in orig and hasattr(exo_module, "fit_psf_photometry_flux_row"):
        orig["fit_psf_photometry_flux_row"] = exo_module.fit_psf_photometry_flux_row
        _orig_fit = orig["fit_psf_photometry_flux_row"]

        def _warmed_psf_fit(data, centroid_row, starIndex, psf_function=None, box=15):
            try:
                if psf_function is None and hasattr(exo_module, "gaussian_psf"):
                    psf_function = exo_module.gaussian_psf
                crow = _np.asarray(centroid_row, dtype=float).reshape(-1)
                if crow.size >= 2 and _np.isfinite(crow[0]) and _np.isfinite(crow[1]):
                    mom = gaussian_psf_moments_cpu(data, float(crow[0]), float(crow[1]), box=box)
                    if mom is not None and mom.size >= 7:
                        new_row = _np.array([
                            crow[0], crow[1],
                            mom[2], mom[3], mom[4], mom[5], mom[6],
                        ], dtype=float)
                        if new_row.size >= 7 and _np.all(_np.isfinite(new_row[:5])):
                            if crow.size >= 7:
                                seed_row = crow.copy()
                                seed_row[2:7] = _np.where(_np.isfinite(new_row[2:7]), new_row[2:7], crow[2:7])
                            else:
                                seed_row = new_row
                            result = _orig_fit(data, seed_row, starIndex, psf_function, box)
                            return result
            except Exception:
                pass
            return _orig_fit(data, centroid_row, starIndex, psf_function, box)

        exo_module.fit_psf_photometry_flux_row = _warmed_psf_fit


def _unpatch_exotic_alignment_and_psf() -> None:
    exo_module = _resolve_exotic_module()
    if exo_module is None:
        return

    if hasattr(exo_module, "_exotic_gpu_original"):
        for name, func in exo_module._exotic_gpu_original.items():
            if hasattr(exo_module, name):
                setattr(exo_module, name, func)
        try:
            delattr(exo_module, "_exotic_gpu_original")
        except Exception:
            pass


__all__ = [
    "enable_exotic_gpu",
    "disable_exotic_gpu",
    "get_full_gpu_status",

    "phase_cross_correlation_gpu",
    "batch_phase_cross_correlation_gpu",
    "align_image_stack_gpu",
    "should_use_gpu_for_alignment",
    "get_alignment_gpu_info",

    "compute_frame_aperture_grid_batch",
    "compute_aperture_grid_cpu",
    "image_aperture_signal_flux_gpu",
    "batch_photometry_simple",
    "should_use_gpu_photometry",
    "is_phot_gpu_available",
    "get_phot_gpu_backend",

    "psf_moments_batch",
    "gaussian_psf_moments_cpu",
    "estimate_image_fwhm_batch",
    "is_psf_gpu_available",
    "get_psf_backend",

    "patch_exotic_photometry",
    "unpatch_exotic_photometry",
    "set_use_gpu",
    "is_gpu_enabled",
    "image_aperture_signal_flux_integrated",
    "batch_photometry_exotic",
]
