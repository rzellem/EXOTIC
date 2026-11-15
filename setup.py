"""Setuptools entry point; package metadata lives in setup.cfg and pyproject.toml.

Build distributions with ``python -m build``. Versioning is generated from Git
by setuptools_scm; publishing is handled by the Python Package workflow.
"""
from setuptools import setup

if __name__ == "__main__":
    setup()
