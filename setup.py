#! /usr/bin/env python3
# ##############################################################################
# SETUPTOOLS PLACEHOLDER
# NOTE: This file supports editable-mode installs using pip3, e.g. 'pip install -e'.
#       It also support building zip-only based source distributions.
# ############################### INSTALL TOOLS  ###############################
# TWINE DOCS: https://twine.readthedocs.io/en/latest/
# pip3 install --upgrade build setuptools_scm twine
# ############################### CHECK VERSION  ###############################
# python3 -m setuptools_scm
# ######################## INSTALL PRODUCT REQUIREMENTS ########################
# pip3 --exists-action w install -r requirements.txt
# ################################### CLEAN  ###################################
# rm -r build dist __pycache__ *.egg* .egg* ; git checkout exotic/version.py ; pip3 uninstall exotic -y
# ############################## RELEASE PRODUCT  ##############################
# (1) manually update exotic/version.py, commit and push
# git add exotic/version.py && git commit -m "#<ticket>: Updated version for release" && git push
# (2) tag using either web ui on github master branch or git cli (here)
# git tag -a -m "#<ticket>: Release version <version>" <version>
# (3) build and package
# git checkout exotic/version.py && python3 -m build --wheel
# git checkout exotic/version.py && python3 setup.py sdist --format=zip
# (3) release and upload to PyPi
# twine check dist/* && twine upload --verbose  dist/*.whl dist/*.zip
# ############################### LOCAL TESTING  ###############################
# python3 -m build --wheel && python3 setup.py sdist --format=zip  # build (zip for source dist)
# pip3 install exotic --no-index --find-links file:///proj/survey-ws/source/EXOTIC/dist/  # install locally


from setuptools import setup

if __name__ == "__main__":
    setup()
