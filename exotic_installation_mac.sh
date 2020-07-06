#!/usr/bin/env sh

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

brew install python

pip3 install --upgrade keyrings.alt
pip3 install -r requirements.txt

chmod 755 exotic.py && python3 exotic.py || python exotic.py
