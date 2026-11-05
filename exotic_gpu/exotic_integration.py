"""
EXOTIC GPU Integration Module
========================
GPU integration helpers for EXOTIC photometry pipeline.

Designed to be **portable across every OS / distribution:
  - No hard-coded user paths (e.g. /home, "C:\\", …)
  - No distribution-specific dependencies (CachyOS, Ubuntu, Windows, macOS, …)
  - Canonical submodule resolution via a single source of truth
    (``EXOTIC.exotic_gpu._module_resolver``).

Canonical / preferred usage (new code in INIT__.py):
    from EXOTIC.exotic_gpu import enable_exotic_gpu
    enable_exotic_gpu()

Legacy usage (kept for backward compatibility):
    from EXOTIC.exotic_gpu.exotic_integration import patch_exotic_photometry
    patch_exotic_photometry()

Author: EXOTIC GPU Optimization Team
Date: 2026-09-01
"""

import numpy as _np
from typing import Tuple, List, Optional, Dict, Any
import warnings as _warnings
import time as _time

# =========================================================
# Cross-platform EXOTIC submodule resolution (DRY).
# Single source of truth: EXOTIC.exotic_gpu._module_resolver.
# =========================================================
try:
    from ._module_resolver import (
        resolve_exotic_module as _resolve_exotic_module_integ,
        EXOTIC_MODULE_CANDIDATES as _EXOTIC_MODULE_CANDIDATES,
    )
except Exception:
    try:
        from EXOTIC.exotic_gpu._module_resolver import (
            resolve_exotic_module as _resolve_exotic_module_integ,
        )
    except Exception:
        # Final fallback: for the user imported exotic_integration.py outside of a full
        # without the EXOTIC directory tree (edge case).
        def _resolve_exotic_module_integ():
            import sys as _sl
            for _cn in ("EXOTIC.exotic.exotic", "exotic.exotic", "EXOTIC.exotic", "exotic"):
                if _cn in _sl.modules:
                    return _sl.modules[_cn]
            try:
                from EXOTIC.exotic import exotic as _e
                return _e
            except ImportError:
                pass
            try:
                import exotic.exotic as _e
                return _e
            except ImportError:
                pass
            try:
                import exotic as _e
                return _e
            except ImportError:
                return None


# =============================================================================
# Photometry GPU/CPU backend (uses modern helpers from photometry_gpu)
# =============================================================================
_HAS_PHOT_GPU = False
_PHOT_BACKEND = "cpu"

try:
    from .photometry_gpu import (
        compute_aperture_grid_cpu,
        compute_frame_aperture_grid_batch,
        should_use_gpu_photometry,
        is_phot_gpu_available,
        get_phot_gpu_backend,
    )
    _PHOT_BACKEND = get_phot_gpu_backend()
    _HAS_PHOT_GPU = is_phot_gpu_available()
except Exception:
    try:
        from EXOTIC.exotic_gpu.photometry_gpu import (
            compute_aperture_grid_cpu,
            compute_frame_aperture_grid_batch,
            should_use_gpu_photometry,
            is_phot_gpu_available,
            get_phot_gpu_backend,
        )
        _PHOT_BACKEND = get_phot_gpu_backend()
        _HAS_PHOT_GPU = is_phot_gpu_available()
    except Exception:
        _warnings.warn("EXOTIC.exotic_gpu.photometry_gpu unavailable; CPU fallback used")

# Global control flag
_use_gpu = _HAS_PHOT_GPU


def set_use_gpu(use_gpu: bool):
    """Toggle GPU usage for the photometry backend only."""
    global _use_gpu
    _use_gpu = use_gpu and _HAS_PHOT_GPU


def is_gpu_enabled() -> bool:
    """Return True if the photometry backend is currently using GPU."""
    return _use_gpu


# =============================================================================
# Safe GPU device information helper (no-op when no GPU)
# =============================================================================
def _get_gpu_info_safe() -> Optional[Dict[str, Any]]:
    """Try to return GPU info; return None on any failure."""
    try:
        from .alignment_gpu import get_gpu_info
        info = get_gpu_info()
        if info.get("available"):
            return info
    except Exception:
        pass
    try:
        import torch  # type: ignore
        if torch.cuda.is_available():
            return {
                "backend": "pytorch",
                "device_name": torch.cuda.get_device_name(0),
                "memory_gb": round(torch.cuda.get_device_properties(0).total_memory / 1e9, 2),
            }
    except Exception:
        pass
    return None


# =============================================================================
# Drop-in compatible helpers for EXOTIC (legacy)
# =============================================================================
def image_aperture_signal_flux(
    image_data: _np.ndarray,
    x_pos: float,
    y_pos: float,
    aperture_radius: float,
    sky_background: Optional[float] = None
) -> float:
    """
    Optimised GPU/CPU version of EXOTIC image_aperture_signal_flux.

    Computes flux within a circular aperture centered on the pixel coordinates
    (x_pos, y_pos).
    """
    try:
        ap = _np.asarray([aperture_radius], dtype=float)
        an = _np.asarray([0.0], dtype=float)
        res = compute_aperture_grid_cpu(image_data, float(x_pos), float(y_pos), ap, an)
        if res is not None and len(res) >= 2:
            flux_grid, bg_grid = res[0], res[1]
            if flux_grid.size > 0:
                flux_val = float(flux_grid[0, 0])
                if sky_background is not None and _np.isfinite(sky_background):
                    bg_val = float(bg_grid[0, 0]) if bg_grid.size > 0 else 0.0
                    if not _np.isfinite(bg_val):
                        bg_val = 0.0
                    if not _np.isfinite(flux_val):
                        flux_val = _image_aperture_signal_flux_cpu_fallback(
                            image_data, x_pos, y_pos, aperture_radius, sky_background
                        )
                    else:
                        # Override the auto background value with the explicit
                        # user-provided one.
                        from math import pi
                        area = pi * aperture_radius * aperture_radius
                        flux_val = flux_val - float(bg_val) * area + sky_background * area
                if not _np.isfinite(flux_val):
                    flux_val = _image_aperture_signal_flux_cpu_fallback(
                        image_data, x_pos, y_pos, aperture_radius, sky_background
                    )
                return flux_val
    except Exception:
        pass
    return _image_aperture_signal_flux_cpu_fallback(
        image_data, x_pos, y_pos, aperture_radius, sky_background
    )


def _image_aperture_signal_flux_cpu_fallback(
    image_data: _np.ndarray,
    x_pos: float,
    y_pos: float,
    aperture_radius: float,
    sky_background: Optional[float] = None
) -> float:
    """Trivial CPU aperture photometry as the final fallback."""
    try:
        height, width = image_data.shape
        r = float(aperture_radius)
        x0 = max(0, int(_np.floor(x_pos - r)))
        x1 = min(width, int(_np.ceil(x_pos + r)) + 1)
        y0 = max(0, int(_np.floor(y_pos - r)))
        y1 = min(height, int(_np.ceil(y_pos + r)) + 1)
        if x0 >= x1 or y0 >= y1:
            return float("nan")
        cutout = _np.asarray(image_data[y0:y1, x0:x1], dtype=float)
        yy, xx = _np.ogrid[y0:y1, x0:x1]
        r2 = (xx - x_pos) ** 2 + (yy - y_pos) ** 2
        mask = r2 <= r * r
        finite = _np.isfinite(cutout)
        valid = mask & finite
        if not _np.any(valid):
            return float("nan")
        flux = float(_np.sum(_np.where(valid, cutout, 0.0)))
        if sky_background is not None and _np.isfinite(sky_background):
            area = float(_np.count_nonzero(valid))
            flux -= float(sky_background) * area
        return flux
    except Exception:
        return float("nan")


def batch_photometry_exotic(
    images: _np.ndarray,
    star_positions: List[Tuple[float, float]],
    aperture_radius: float,
    sky_backgrounds: Optional[_np.ndarray] = None
) -> _np.ndarray:
    """
    Batch aperture photometry across multiple images and multiple stars.

    Args:
        images: (N, H, W) or (H, W) frames to process
        star_positions: list of (x, y) pixel coordinates per star
        aperture_radius: radius in pixels (same for all apertures)
        sky_backgrounds: optional (N,) per-frame sky background values

    Returns:
        fluxes: (N, n_stars) numpy array of measured fluxes
    """
    images = _np.asarray(images)
    if images.ndim == 2:
        images = images[_np.newaxis, :, :]
    n_images = images.shape[0]
    n_stars = len(star_positions)
    fluxes = _np.full((n_images, n_stars), _np.nan, dtype=float)
    ap = _np.asarray([aperture_radius], dtype=float)
    an = _np.asarray([0.0], dtype=float)

    try:
        coords = [(float(x), float(y)) for (x, y) in star_positions]
        use_batch_gpu = _use_gpu and _HAS_PHOT_GPU and should_use_gpu_photometry(n_stars * 1, pixels=images[0].size if n_images else 1)
    except Exception:
        use_batch_gpu = False

    for i in range(n_images):
        frame = images[i]
        sky_i = float(sky_backgrounds[i]) if (sky_backgrounds is not None and i < len(sky_backgrounds) and _np.isfinite(sky_backgrounds[i])) else None
        if use_batch_gpu and i == 0:
            try:
                batch_res = compute_frame_aperture_grid_batch(frame, coords, ap, an, use_gpu_if_available=True)
                for j in range(n_stars):
                    try:
                        fg = batch_res[j][0]
                        if fg.size > 0 and _np.isfinite(fg[0, 0]):
                            fv = float(fg[0, 0])
                            if sky_i is not None:
                                from math import pi
                                area = pi * aperture_radius * aperture_radius
                                fv = fv + sky_i * area
                            fluxes[i, j] = fv
                        else:
                            fluxes[i, j] = image_aperture_signal_flux(frame, coords[j][0], coords[j][1], aperture_radius, sky_i)
                    except Exception:
                        fluxes[i, j] = image_aperture_signal_flux(frame, coords[j][0], coords[j][1], aperture_radius, sky_i)
                continue
            except Exception:
                pass
        for j in range(n_stars):
            x, y = star_positions[j]
            fluxes[i, j] = image_aperture_signal_flux(frame, x, y, aperture_radius, sky_i)

    return fluxes


# =============================================================================
# EXOTIC patching helpers (legacy; enable_exotic_gpu is the preferred entry point)
# =============================================================================
def patch_exotic_photometry():
    """
    Monkey-patch EXOTIC photometry functions with optimised variants.

    NOTE: The preferred and most complete entry point is:
        from EXOTIC.exotic_gpu import enable_exotic_gpu
        enable_exotic_gpu()
    """
    exo_module = _resolve_exotic_module_integ()
    if exo_module is None:
        _warnings.warn("Could not locate the EXOTIC module to patch photometry.")
        return

    try:
        if not hasattr(exo_module, "_exotic_gpu_integ_orig"):
            exo_module._exotic_gpu_integ_orig = {}
        orig = exo_module._exotic_gpu_integ_orig

        if hasattr(exo_module, "image_aperture_signal_flux") and "image_aperture_signal_flux" not in orig:
            orig["image_aperture_signal_flux"] = exo_module.image_aperture_signal_flux
            exo_module.image_aperture_signal_flux = image_aperture_signal_flux

        status = "GPU" if (_use_gpu and _HAS_PHOT_GPU) else "CPU (optimised)"
        print(f"OK  EXOTIC photometry patched: {status}")
        if _use_gpu and _HAS_PHOT_GPU:
            info = _get_gpu_info_safe()
            if info is not None:
                print(f"     Backend: {info.get('backend', '?').upper()}  Device: {info.get('device_name', 'GPU')}")
                if "memory_gb" in info:
                    print(f"     VRAM: {info['memory_gb']} GB")
    except Exception as e:
        _warnings.warn(f"Error patching EXOTIC photometry: {e}")


def unpatch_exotic_photometry():
    """Restore EXOTIC's original photometry functions (if patched previously)."""
    exo_module = _resolve_exotic_module_integ()
    if exo_module is None:
        return
    try:
        if hasattr(exo_module, "_exotic_gpu_integ_orig"):
            orig = exo_module._exotic_gpu_integ_orig
            for name, func in orig.items():
                if hasattr(exo_module, name):
                    setattr(exo_module, name, func)
            try:
                delattr(exo_module, "_exotic_gpu_integ_orig")
            except Exception:
                pass
            print("OK  EXOTIC photometry restored to original")
    except Exception:
        pass


# =============================================================================
# Legacy entry point for testing
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("EXOTIC GPU Integration Module — Status")
    print("=" * 70)
    print(f"  Photometry backend: {_PHOT_BACKEND.upper()}"
          + (" (GPU)" if _HAS_PHOT_GPU else " (CPU fallback)"))
    print(f"  GPU enabled (flag): {_use_gpu}")
    info = _get_gpu_info_safe()
    if info is not None:
        print(f"  GPU device:         {info.get('device_name', '?')}")
        if "memory_gb" in info:
            print(f"  VRAM available:     {info['memory_gb']} GB")
    print()
    print("  (Recommended usage: from EXOTIC.exotic_gpu import enable_exotic_gpu)")
    print("=" * 70)
