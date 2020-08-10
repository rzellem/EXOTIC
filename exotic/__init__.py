import importlib_metadata as metadata

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

