#!/bin/sh

git clone https://github.com/rzellem/EXOTIC.git

sudo apt-get -y update && sudo apt-get -y upgrade
apt-get -y install python3-pip

pip3 install numpy
pip3 install astropy
pip3 install python-dateutil
pip3 install barycorrpy
pip3 install --upgrade keyrings.alt
pip3 install matplotlib
pip3 install pymc3
pip3 install theano
pip3 install photutils
apt-get -y install g++
python3 ./EXOTIC/exotic.py
