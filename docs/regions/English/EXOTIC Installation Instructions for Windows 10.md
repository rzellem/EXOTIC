# Install EXOTIC on Windows

EXOTIC requires **Python 3.12 or newer**. The repository includes Windows x64
UltraNest wheels for CPython 3.12, 3.13, and 3.14, plus a SEP wheel for 3.14.
The instructions below use the latest Python 3.14 maintenance release.
See the [wheelhouse README](../../../wheelhouse/README.md) for details.

## Create a Python environment

Install the [Python install manager](https://www.python.org/downloads/windows/)
if you do not already have it. Open **Command Prompt** and install or update the
standard x64 Python 3.14 runtime:

```bat
pymanager install --update 3.14-64
```

This selects the latest 3.14 maintenance release rather than whichever Python
major/minor release is newest. `pymanager` selects the install manager even when
the older `py` launcher is also installed. You can alternatively install the
latest 3.14 **Windows installer (64-bit)** from
[Python's release list](https://www.python.org/downloads/windows/).

Create a dedicated environment and activate it:

```bat
pymanager exec -V:3.14-64 -m venv "%USERPROFILE%\.venvs\exotic314"
"%USERPROFILE%\.venvs\exotic314\Scripts\activate.bat"
python --version
python -m pip install --upgrade pip
```

Activating the environment makes its Python and installed commands available in
that terminal. If you used the traditional installer and older launcher, create
the environment with `py -3.14-64 -m venv "%USERPROFILE%\.venvs\exotic314"`
instead. For Python 3.12 or 3.13, replace `3.14-64` with `3.12-64` or `3.13-64`
and use a separate environment directory. The bundled wheels target standard
CPython on x64; ARM64 and free-threaded Python (`3.14t`) use different wheel tags.

## Install and launch EXOTIC

Download and extract the [EXOTIC source repository](https://github.com/rzellem/EXOTIC/archive/refs/heads/main.zip),
or use an existing checkout. Keep the `wheelhouse` folder alongside the launcher.
In the same Command Prompt, change to the extracted repository directory. To
install this checkout, including its current Python 3.14 support, run:

```bat
python -m pip install --find-links wheelhouse --only-binary=ultranest -e .
python -m pip check
exotic --help
exotic-gui
```

The pip command offers the bundled wheels and requires a binary UltraNest
package, avoiding a local UltraNest source build. The bundled SEP wheel also
avoids compiling that dependency on Python 3.14. `exotic` is the command-line
interface; `exotic-gui` opens the graphical interface supplied by the same
environment.

Alternatively, run `run_exotic_windows.bat` to install or update the **published
EXOTIC package** and then open its GUI. The launcher uses the active environment's
Python for both pip and the GUI. Installation errors stop the launch and are
shown in the terminal.

## Run EXOTIC again

In a new Command Prompt, activate the same environment:

```bat
"%USERPROFILE%\.venvs\exotic314\Scripts\activate.bat"
exotic-gui
```

Use `exotic` for the command-line interface. For an editable checkout install,
source edits are used directly; rerun its pip install command when dependencies
change. The batch launcher updates to the published package instead.

The activation commands above are for Command Prompt. In PowerShell, activate
with `& "$env:USERPROFILE\.venvs\exotic314\Scripts\Activate.ps1"`. You can also
use the environment without activation or changing PowerShell execution policy:

```powershell
& "$env:USERPROFILE\.venvs\exotic314\Scripts\python.exe" -m pip install --find-links wheelhouse --only-binary=ultranest -e .
& "$env:USERPROFILE\.venvs\exotic314\Scripts\exotic-gui.exe"
```

The [repository README](../../../README.md) also covers Conda, Google Colab,
sample data, and optional MPI support.

See [Python's Windows guide](https://docs.python.org/3.14/using/windows.html)
for install manager, runtime selection, and virtual environment details.
