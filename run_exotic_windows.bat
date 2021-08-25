ECHO OFF

pip install wheel
pip install setuptools
pip install --upgrade python-certifi-win32
pip install --upgrade exotic
start cmd /k exotic-gui
