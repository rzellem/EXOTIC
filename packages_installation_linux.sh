#!/bin/sh

sudo apt update
sudo apt -y upgrade
sudo apt -y install python3-pip

pip3 install numpy
pip3 install astropy
sudo pip3 install python-dateutil
pip3 install barycorrpy
pip3 install --upgrade keyrings.alt
pip3 install matplotlib
pip3 install pymc3
pip3 install theano
pip3 install photutils
sudo apt-get -y install g++
python3 exotic.py
