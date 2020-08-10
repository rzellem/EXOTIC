import codecs
from pathlib import Path


def read_relative(rel_path):
    here = Path(__file__).parent
    with codecs.open(here / rel_path, 'r') as fp:
        return fp.read()


def version_read(rel_path):
    version_semantic_fallback = "0.1.0"
    for line in read_relative(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError(f"ERROR: Unable to parse version string from {rel_path}. "
                           f"Using fallback of {version_semantic_fallback}...")
