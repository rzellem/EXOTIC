@ECHO OFF
SETLOCAL

REM Use the active environment and run this checkout directly.
python -c "import sys; raise SystemExit(sys.version_info < (3, 12))"
if errorlevel 1 (
    echo ERROR: Python 3.12 or newer is required. Activate your EXOTIC environment or install a supported Python on PATH.
    exit /b 65
)
cd /d "%~dp0"
if errorlevel 1 exit /b %errorlevel%
python -m exotic.api.native_build
if errorlevel 1 exit /b %errorlevel%
python -c "from exotic.exotic_gui import main; main()"
exit /b %errorlevel%
