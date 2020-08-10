import importlib_metadata as metadata

from .api.version import version_read

ignore = True

__version__ = "0.1.0"  # placeholder for error in reporting true version
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

