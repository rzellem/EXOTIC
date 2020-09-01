#! /usr/bin/env python3
# python3 setup.py --version
# FOR RELEASE:
# manually update exotic/version.py
# git add exotic/version.py && git commit -m "#<ticket>:Updated version for release" && git push
# release/tag
# git checkout exotic/version.py && python3 setup.py sdist bdist_wheel --universal  # upload
# pip3 install twine  # https://twine.readthedocs.io/en/latest/
# twine check dist/* && twine upload dist/*
# LOCAL TESTING:
# pip3 install --exists-action w --progress-bar ascii --use-feature=2020-resolver -r requirements.txt
# pip3 install exotic --no-index --find-links file:///proj/survey-ws/source/EXOTIC/dist/
# rm -r dist && pip3 uninstall exotic

import re
import time
from pathlib import Path

import setuptools

# Package meta-data sane defaults
AUTHOR = "Exoplanet Watch at NASA JPL"
AUTHOR_EMAIL = "exoplanetwatch@jpl.nasa.gov"
DESCRIPTION = "EXOTIC: EXOplanet Transit Interpretation Code"
NAME = "exotic"
PYTHON_REQUIREMENTS = ">=3.7.0"
URL = "https://github.com/rzellem/EXOTIC"

REQUIREMENTS_SETUP = ['setuptools_scm',
                      'importlib-metadata ~= 1.7 ; python_version >= "3.7"']


def description_read():
    description_long = ""
    description_path = Path("README.md")
    if not description_path.exists():
        return description_long
    with description_path.open('r') as f:
        description_long = f.read()
    return description_long


def license_read():
    annum = time.gmtime().tm_year
    # provide one-line summary per LICENSE spec
    lic = f"Proprietary -- Copyright (c) 2019-{annum}, California Institute of Technology." 
    return lic


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


setuptools.setup(name=NAME,
                 use_scm_version={
                     'write_to': 'exotic/version.py',
                     'write_to_template': '__version__ = "{version}"'
                 },
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
                     'Intended Audience :: End Users/Desktop',
                     'Intended Audience :: Science/Research',
                     'License :: Other/Proprietary License',
                     'Programming Language :: Python',
                     'Programming Language :: Python :: 3',
                     'Programming Language :: Python :: 3.7',
                     'Programming Language :: Python :: Implementation :: CPython',
                     'Programming Language :: Python :: Implementation :: PyPy',
                     'Topic :: Scientific/Engineering :: Astronomy'
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
                 setup_requires=REQUIREMENTS_SETUP,
                 entry_points={
                     'console_scripts': [
                         'exotic = exotic.exotic:main',
                     ],
                 }
                 )
