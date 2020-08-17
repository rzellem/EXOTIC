import importlib_metadata as metadata
from pathlib import Path
import sys

# extend PYTHONPATH
path_update = ['.', str(Path(__file__).resolve().parent)]
for p in path_update:
    if p not in sys.path:
        sys.path.append(p)

try:  # module import
    from .api.versioning import version_read
except ImportError:  # package import
    from api.versioning import version_read

ignore = True

try:
    __version__ = metadata.version(__name__)
except metadata.PackageNotFoundError:
    # package is not installed, try reading from exotic script
    try:
        __version__ = version_read("exotic.py")
    except IOError:
        # unable to read from exotic script
        __version__ = "unknown"
        pass
    pass
