@ECHO OFF
SETLOCAL

set "WHEELHOUSE=%~dp0wheelhouse"
if not exist "%WHEELHOUSE%" (
    echo ERROR: The local UltraNest wheelhouse was not found: "%WHEELHOUSE%"
    exit /b 1
)

REM Use Python from the active environment for pip and the GUI.
python -c "import sys; raise SystemExit(sys.version_info < (3, 12))"
if errorlevel 1 (
    echo ERROR: Python 3.12 or newer is required. Activate your EXOTIC environment or install a supported Python on PATH.
    exit /b 65
)

python -m pip --version
if not errorlevel 1 goto install_exotic
python -m ensurepip --upgrade
if errorlevel 1 exit /b %errorlevel%

:install_exotic
REM Pip manages build dependencies from the package metadata.
python -m pip install --upgrade --find-links "%WHEELHOUSE%" --only-binary ultranest "ultranest>=4.5.0" exotic
if errorlevel 1 exit /b %errorlevel%

REM Use the installed package rather than a checkout in the current directory.
python -P -c "from exotic.exotic_gui import main; main()"
exit /b %errorlevel%
