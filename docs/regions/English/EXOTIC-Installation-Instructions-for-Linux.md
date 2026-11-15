# Install EXOTIC on Linux

EXOTIC requires Python 3.12 or newer. The commands below use Bash and Python
3.14. Install the latest 3.14 maintenance release from your Python provider.
For Python 3.12 or 3.13, substitute `python3.12` or `python3.13`; `python3` also
works when it selects a supported interpreter.

## Python virtual environment

With your selected Python interpreter installed:

```bash
python3.14 --version
python3.14 -m venv ~/.venvs/exotic314
source ~/.venvs/exotic314/bin/activate
python -m pip install --upgrade pip
python -m pip install exotic
```

`venv` creates the environment and installs pip in it. Activating it selects
that environment's Python and installed commands. EXOTIC and its dependencies
are installed in the environment.

To use the source checkout instead of the published package, run
`python -m pip install -e .` from the EXOTIC repository root in the active
environment, in place of `python -m pip install exotic`.

## Conda environment

If you already use Miniconda or Anaconda, create a dedicated environment:

```bash
conda create -n exotic python=3.14 pip
conda activate exotic
python -m pip install exotic
```

This selects the supported Python version explicitly instead of cloning the
base environment's interpreter and packages.

## Run EXOTIC

Activate your environment, then launch the command-line interface:

```bash
exotic
```

For the optional graphical interface, install the Tkinter package matching
your Python interpreter through your distribution's package manager, then run:

```bash
exotic-gui
```

In a new terminal, activate the same environment with
`source ~/.venvs/exotic314/bin/activate` or `conda activate exotic` before running
EXOTIC. Leave it with `deactivate` or `conda deactivate`, respectively.

See the [repository README](../../../README.md) for optional MPI support,
sample data, and reduction instructions.

References: [Python venv documentation](https://docs.python.org/3.14/library/venv.html)
and [Python Packaging User Guide](https://packaging.python.org/en/latest/tutorials/installing-packages/).
