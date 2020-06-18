#!/bin/sh

git clone https://github.com/rzellem/EXOTIC.git

sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get -y update && sudo apt-get -y upgrade
sudo apt install python3
sudo apt-get -y install python3-pip

pip3 install numpy
pip3 install astropy
pip3 install python-dateutil
pip3 install barycorrpy
pip3 install --upgrade keyrings.alt
pip3 install matplotlib
pip3 install pymc3
pip3 install theano
pip3 install photutils
sudo apt-get -y install g++
pip3 install astroalign

python3 ./EXOTIC/exotic.py
