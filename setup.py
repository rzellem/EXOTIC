#! /usr/bin/env python3
# INSTALL TOOLS:
# pip3 install -r requirements.txt ; pip3 install --upgrade build twine  # https://twine.readthedocs.io/en/latest/
# CLEAN PREVIOUS BUILDS:
# rm -r build dist __pycache__ *.egg* .egg* ; pip3 uninstall exotic ; git checkout exotic/version.py # pipenv uninstall exotic
# CHECK VERSION:
# python3 setup.py --version
# RELEASE PRODUCT:
# ############################## RELEASE PRODUCT  ##############################
# (1) manually update exotic/version.py
# git add exotic/version.py && git commit -m "#<ticket>: Updated version for release" && git push
# (2) tag using web ui on github master branch
# or, alternately, git tag -a -m "#<ticket>: Release version <version>" <version>
# (3) build and package
# git checkout exotic/version.py && python3 -m build --sdist
# python3 setup.py sdist --format=gztar,zip bdist_wheel --universal  # build
# (3) release
# twine check dist/* && twine upload dist/*  # upload
# ############################### LOCAL TESTING  ###############################
# python3 setup.py sdist --format=gztar,zip bdist_wheel --universal  # build
# pip3 --exists-action w install -r requirements.txt  # install req's
# pip3 install exotic --no-index --find-links file:///proj/survey-ws/source/EXOTIC/dist/  # install locally


from setuptools import setup

if __name__ == "__main__":
    setup()
