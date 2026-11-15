# Contributing to EXOTIC

EXOTIC is an open source project that welcomes contributions. Please fork the repository and submit a pull request to
the develop branch for your addition(s) to be reviewed.

## Python and Astropy compatibility

The test workflow runs on Windows and Linux with Python 3.12, 3.13, and 3.14.
Each Python version tests both the declared Astropy minimum (8.0.1) and the
latest compatible Astropy release. Each job installs EXOTIC and the selected
Astropy together through pip's resolver, runs `pip check`, checks runtime
imports without MPI, and runs the full test suite. Windows jobs explicitly
install the bundled UltraNest wheel matching their Python ABI, plus the bundled
SEP wheel for Python 3.14. The native dependency tests execute compiled
pair-distance and source-extraction operations as well as checking imports.

To reproduce one combination, activate an environment using the desired Python
version and run from the repository root:

```bash
python -m pip install --upgrade -e ".[test]" "astropy==8.0.1"
python -m pip check
python -m pytest -q
```

On Windows, add `--find-links wheelhouse --only-binary=ultranest` to the install
command. To test the latest compatible Astropy instead of the minimum, replace
`astropy==8.0.1` with `astropy>=8.0.1`. The Windows wheels cover standard x64
CPython 3.12, 3.13, and 3.14; see the [wheelhouse notes](wheelhouse/README.md).

When changing the minimum in `requirements.txt`, update the minimum matrix
selection in `.github/workflows/python-tests.yml` and these reproduction commands
so CI continues to test the actual supported floor.

### Local verification, 2026-10-24

Windows x64, CPython 3.14.9, NumPy 2.5.4, and the bundled UltraNest 4.5.0 and
SEP 1.4.1 builds:

| Astropy | Passing tests |
| --- | --- |
| 8.0.1 | 1,896 full-suite tests, plus the two native dependency tests separately |
| 8.0.2 | 1,898 full-suite tests, including the native dependency tests |

Both installations passed `pip check`. CLI help, CLI/GUI module imports, Tk
initialization, and a small native UltraNest sampling run also passed. A fresh
pip resolution selected both bundled CPython 3.14 wheels. The remote Windows
and Linux matrix runs are pending.
