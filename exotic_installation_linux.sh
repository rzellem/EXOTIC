#!/usr/bin/env sh

sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get -y update && sudo apt-get -y upgrade
sudo apt-get install python3
sudo apt-get -y install python3-pip
sudo apt-get -y install g++
pip3 install --upgrade pip


pip3 install -r requirements.txt
pip3 install --upgrade keyrings.alt

chmod 755 exotic/exotic.py && python3 exotic/exotic.py || python exotic/exotic.py
