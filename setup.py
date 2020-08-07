#! /usr/bin/env python3
# python3 setup.py sdist bdist_wheel upload

from pathlib import Path
import re
import setuptools
import time


# Package meta-data.
VERSION_SEMANTIC_FALLBACK = "0.1.0"
VERSION = VERSION_SEMANTIC_FALLBACK
AUTHOR = "Exoplanet Watch at NASA JPL"
AUTHOR_EMAIL = "jengelke@jpl.nasa.gov"
DESCRIPTION = "EXOTIC: EXOplanet Transit Interpretation Code"
NAME = "EXOTIC"
PYTHON_REQUIREMENTS = ">=3.7.0"
URL = "https://github.com/rzellem/EXOTIC"

REQUIREMENTS_SETUP = ['setuptools_scm',
                      'importlib-metadata ~= 1.0 ; python_version >= "3.7"']


def license_read():
    annum = time.gmtime().tm_year
    # provide basic default in case LICENSE file missing
    lic = f"""
**      Copyright (c) 2019-{annum}, California Institute of Technology. 
**    All rights reserved.  Based on Government Sponsored Research under 
**    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001. 
"""
    lic_path = Path("LICENSE")
    if not lic_path.exists():
        return lic
    with lic_path.open('r') as f:
        lic = f.read()
    lic = re.sub(r"\s*-\s*\d{4}\s*,", f"-{annum},", lic, 0, re.MULTILINE)
    return lic


def description_read():
    description_long = ""
    description_path = Path("README.md")
    if not description_path.exists():
        return description_long
    with description_path.open('r') as f:
        description_long = f.read()
    return description_long


def requirements_read():
    requirements = []
    requirements_path = Path("requirements.txt")
    if not requirements_path.exists():
        return requirements
    with requirements_path.open('r') as f:
        for line in f:
            # detect and skip comment lines
            if re.match(r"^\s*#.*", line):
                continue
            requirements.append(str(line).strip())
    return requirements


# TODO: Add scripts, console_scripts or entry_points to execute exotic
setuptools.setup(name=NAME,
                 use_scm_version=True,
                 description=DESCRIPTION,
                 long_description=description_read(),
                 long_description_content_type='text/markdown',
                 url=URL,
                 author=AUTHOR,
                 author_email=AUTHOR_EMAIL,
                 license=license_read(),
                 classifiers=[
                     # Trove classifiers
                     # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
                     'Development Status :: 4 - Beta',
                     'License :: Other/Proprietary License',
                     'Programming Language :: Python',
                     'Programming Language :: Python :: 3',
                     'Programming Language :: Python :: 3.7',
                     'Programming Language :: Python :: Implementation :: CPython',
                     'Programming Language :: Python :: Implementation :: PyPy',
                     'Scientific/Engineering :: Astronomy'
                 ],
                 keywords='nasa jpl exoplanet transit citizen science astronomy bayesian nested-sampler',
                 project_urls={
                    "Documentation": "https://github.com/rzellem/EXOTIC/wiki",
                    "Site": "https://exoplanets.nasa.gov/exoplanet-watch",
                    "Source": "https://github.com/rzellem/EXOTIC",
                    "Tracker": "https://github.com/rzellem/EXOTIC/issues"
                 },
                 #  https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
                 packages=setuptools.find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
                 include_package_data=True,
                 zip_safe=False,
                 install_requires=requirements_read(),
                 python_requires=PYTHON_REQUIREMENTS,
                 setup_requires=REQUIREMENTS_SETUP
                 )
