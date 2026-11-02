"""
EXOTIC GPU — Cross-platform module resolver
============================================
Shared utility used to locate the EXOTIC main package and its submodules
independently of how the user invokes the code:
  - installed package:     `from EXOTIC.exotic.exotic import main`
  - direct script launch:  `python EXOTIC/exotic/exotic.py`
  - local dev tree:        `cd EXOTIC/.. ; python -m EXOTIC.exotic.exotic`
  - any operating system (Linux / Windows / macOS).

There are zero hard-coded user paths, distribution-specific calls or
GPU-model assumptions here. Anything platform-specific is wrapped in a
try/except so it can degrade gracefully.
"""

from __future__ import annotations

import sys as _sys
import importlib.util as _ilu
from typing import Optional, Any


# Canonical list of module names to try, ordered by preference. The first
# entry that exists (either loaded in sys.modules or importable) wins.
#
# Target repository layout (official EXOTIC repo):
#   EXOTIC/            (repo root)
#     ├── exotic/      (main package, untouched by GPU rev)
#     │   └── exotic.py
#     └── exotic_gpu/  (NEW GPU subpackage — sibling of exotic/)
EXOTIC_MODULE_CANDIDATES = (
    "EXOTIC.exotic.exotic",    # Canonical shape of the official repository.
    "exotic.exotic",           # User imports it as top-level "exotic".
    "EXOTIC.exotic",           # Rare install without the exotic.exotic alias.
    "exotic",                  # Minimal script-only usage.
)


def _iter_sys_modules(keys):
    """Safe copy of a module-key iterable (avoids RuntimeError if mutated)."""
    return list(keys)


def resolve_exotic_module() -> Optional[Any]:
    """Locate the ``exotic.exotic`` submodule — the main EXOTIC module that
    contains the public functions we monkey-patch (``phase_cross_correlation``,
    ``compute_star_aperture_grid``, ``fit_psf_photometry_flux_row``, …).

    Returns the module object, or ``None`` if it cannot be located through
    any mechanism. Never raises; internal errors are caught and converted
    to a ``None`` return so the caller can fall back.

    Algorithm (in priority order):
      1. Search :data:`sys.modules` first (candidates that were already
         imported earlier). Important because ``EXOTIC.exotic`` often loads
         partially before the monkey-patch is applied.
      2. ``importlib.import_module`` for each candidate name.
      3. Filesystem search if the current working directory contains an
         ``EXOTIC/exotic/`` or ``exotic/`` tree (local dev mode without a
         prior ``pip install -e``). Only used when steps 1 and 2 both fail.
    """
    # Step 1: already present in sys.modules
    try:
        modules = _iter_sys_modules(_sys.modules.keys())
        for name in EXOTIC_MODULE_CANDIDATES:
            if name in modules:
                m = _sys.modules[name]
                if m is not None and hasattr(m, "__name__"):
                    return m
    except Exception:
        pass

    # Step 2: direct import by canonical name
    import importlib
    for name in EXOTIC_MODULE_CANDIDATES:
        try:
            return importlib.import_module(name)
        except ImportError:
            continue
        except Exception:
            # Non-ImportError failures (circular import, etc.) do not
            # block the remaining candidates.
            continue

    # Step 3: local development tree without pip install — look for
    # "EXOTIC/exotic/exotic.py" or "exotic/exotic.py" on sys.path entries.
    import os.path as _op
    for base in list(_sys.path) + [_sys.path[0] if _sys.path else "", ""]:
        try:
            if not base:
                import os
                base = os.getcwd()
            for sub in ("EXOTIC/exotic", "exotic"):
                py_file = _op.join(base, sub, "exotic.py")
                if _op.isfile(py_file) and _op.isdir(_op.join(base, sub)):
                    try:
                        # Use spec_from_file_location to avoid name collisions
                        spec = _ilu.spec_from_file_location(
                            sub.replace("/", ".") + ".exotic",
                            py_file,
                            submodule_search_locations=[_op.join(base, sub)],
                        )
                        if spec is None or spec.loader is None:
                            continue
                        mod = _ilu.module_from_spec(spec)
                        key1 = "EXOTIC.exotic.exotic" if sub.startswith("EXOTIC") else "exotic.exotic"
                        _sys.modules.setdefault(key1, mod)
                        spec.loader.exec_module(mod)
                        return mod
                    except Exception:
                        continue
        except Exception:
            continue

    return None


def resolve_exotic_exotic_gpu_submodule(submodule: str) -> Optional[Any]:
    """Resolve an ``EXOTIC.exotic_gpu`` submodule (``alignment_gpu``,
    ``photometry_gpu``, ``psf_fitting_gpu``, ``exotic_integration``, …).

    Equivalent to ``import EXOTIC.exotic_gpu.<submodule>`` but with a
    fallback to ``import exotic_gpu.<submodule>`` for users who have the
    tree exposed as a top-level package instead of the canonical
    ``EXOTIC.exotic_gpu`` shape.
    """
    import importlib
    for prefix in ("EXOTIC.exotic_gpu", "exotic_gpu"):
        try:
            return importlib.import_module(f"{prefix}.{submodule}")
        except ImportError:
            continue
        except Exception:
            continue
    return None


__all__ = [
    "EXOTIC_MODULE_CANDIDATES",
    "resolve_exotic_module",
    "resolve_exotic_exotic_gpu_submodule",
]
