#!/bin/sh

git clone https://github.com/rzellem/EXOTIC.git

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
##export PATH="/usr/local/opt/python/libexec/bin:$PATH"

brew install python
brew install pip

pip install numpy
pip install astropy
pip install python-dateutil
pip install barycorrpy
pip install --upgrade keyrings.alt
pip install matplotlib
pip install pymc3
pip install theano
pip install photutils

#pip -y install g++
python ./EXOTIC/exotic.py
