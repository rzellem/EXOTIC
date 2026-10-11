@ECHO OFF

set "WHEELHOUSE=%~dp0wheelhouse"
if not exist "%WHEELHOUSE%" (
    echo ERROR: The local UltraNest wheelhouse was not found: "%WHEELHOUSE%"
    exit /b 1
)

pip --version
pip install "importlib-metadata>=3.6;python_version<='3.7'" oldest-supported-numpy "setuptools>=62.6" "setuptools_scm[toml]>=6.4.2" "wheel>=0.37.1"
REM Optional on older Windows certificate setups:
REM pip install --upgrade python-certifi-win32
pip install --upgrade --find-links "%WHEELHOUSE%" --only-binary ultranest "ultranest>=4.5.0" exotic

start cmd /k exotic-gui
