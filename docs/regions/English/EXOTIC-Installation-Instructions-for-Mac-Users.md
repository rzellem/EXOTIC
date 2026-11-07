# Install EXOTIC on macOS

EXOTIC requires **Python 3.12 or newer**. The commands below use Python 3.14.
Use an existing supported Python or install the latest 3.14 maintenance release from
[Python for macOS](https://www.python.org/downloads/macos/).

## Create an environment and install EXOTIC

Open Terminal and run:

```bash
python3.14 -m venv ~/.venvs/exotic314
source ~/.venvs/exotic314/bin/activate
python -m pip install --upgrade pip
python -m pip install exotic
```

For Python 3.12 or 3.13, substitute `python3.12` or `python3.13`; `python3` also
works when it selects a supported interpreter. The environment keeps EXOTIC's
dependencies together and ensures pip uses the same Python as EXOTIC.

To use a source checkout instead of the published package, run
`python -m pip install -e .` from the repository root in the active environment,
in place of `python -m pip install exotic`.

## Run EXOTIC

For the graphical interface:

```bash
exotic-gui
```

Use `exotic` for the command-line interface. The GUI uses Tkinter; if your Python
provider packages Tk separately, install its matching Tk support for the GUI.

In a new Terminal window, activate the environment again with
`source ~/.venvs/exotic314/bin/activate`. Leave it with `deactivate`.

If you use a repository checkout, you can also run `./run_exotic_mac.command`
from its root with the environment active. This installs or updates the
published EXOTIC package before opening the GUI. Installation errors stop the
launch and remain visible in the terminal.

See the [repository README](../../../README.md) for Conda, Google Colab
tutorials, sample data, and optional MPI support.
