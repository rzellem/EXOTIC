import importlib_metadata as metadata

ignore = True

__version__ = "0.1.0"  # placeholder for error in reporting true version
try:
    __version__ = metadata.version(__name__)
except metadata.PackageNotFoundError:
    # package is not installed
    pass

