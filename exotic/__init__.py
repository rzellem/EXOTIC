# ########################################################################### #
#    Copyright (c) 2019-2020, California Institute of Technology.
#    All rights reserved.  Based on Government Sponsored Research under
#    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:
#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#      3. Neither the name of the California Institute of
#         Technology (Caltech), its operating division the Jet Propulsion
#         Laboratory (JPL), the National Aeronautics and Space
#         Administration (NASA), nor the names of its contributors may be
#         used to endorse or promote products derived from this software
#         without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
#    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ########################################################################### #
#    EXOplanet Transit Interpretation Code (EXOTIC)
#    # NOTE: See companion file version.py for version info.
# ########################################################################### #

# Import metadata handling based on Python version
try:
    from importlib import metadata  # Python 3.8+
except ImportError:
    # Fallback to importlib_metadata for older versions
    import importlib_metadata as metadata  # Python <3.8

from pathlib import Path
from importlib import import_module
import sys

# Extend PYTHONPATH to include current directory and parent directory
# This ensures imports can find modules in the current package structure
path_update = ['.', str(Path(__file__).resolve().parent)]
for p in path_update:
    if p not in sys.path:
        sys.path.append(p)

# Import custom version reading function for fallback scenario
try:  
    # First attempt: import as a relative module (when running as part of a package)
    from .api.versioning import version_read
except ImportError:
    # Second attempt: import as an absolute package (when running standalone)
    from api.versioning import version_read

# Version detection 
try:
    # First attempt: get version from package metadata
    __version__ = metadata.version(__name__)
except metadata.PackageNotFoundError:
    # Second attempt: package is not installed, try reading version from exotic script
    try:
        __version__ = version_read("exotic.py")
    except IOError:
        # Unable to read from exotic script
        __version__ = "unknown"


def _load_runtime_callable(name):
    """Load CLI entry points lazily without importing the full runtime at package import."""
    current_module = sys.modules.get(__name__)
    module_names = [f"{__name__}.exotic", "exotic.exotic"]

    seen = set()
    for module_name in module_names:
        if module_name in seen:
            continue
        seen.add(module_name)

        try:
            module = import_module(module_name)
        except ModuleNotFoundError as exc:
            if exc.name == module_name:
                continue
            raise
        if module is current_module:
            continue
        runtime_callable = getattr(module, name, None)
        if runtime_callable is not None:
            return runtime_callable

    raise ImportError(f"cannot import name '{name}' from '{__name__}'")


def main(*args, **kwargs):
    return _load_runtime_callable("main")(*args, **kwargs)


def cli(*args, **kwargs):
    return _load_runtime_callable("cli")(*args, **kwargs)


__all__ = ["__version__", "main", "cli"]
