ECHO OFF

pip --version
pip install "importlib-metadata>=3.6;python_version<='3.7'" oldest-supported-numpy setuptools>=62.6 setuptools_scm[toml]>=6.4.2 wheel>=0.37.1
# pip install --upgrade python-certifi-win32
pip install --upgrade exotic

start cmd /k exotic-gui
