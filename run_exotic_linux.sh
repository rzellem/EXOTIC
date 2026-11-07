#!/usr/bin/env bash

# Use the same interpreter for installation and the GUI, including an active
# virtual environment or Conda environment. The minimum matches setup.cfg.
py_runner=""
for app in python3 python; do
    if command -v "$app" >/dev/null 2>&1 && \
        "$app" -c 'import sys; raise SystemExit(sys.version_info < (3, 12))'; then
        py_runner="$app"
        break
    fi
done

if [[ -z "$py_runner" ]]; then
    echo "ERROR: Python 3.12 or newer is required. Activate your EXOTIC environment or install a supported Python."
    exit 65
fi

"$py_runner" --version
if ! "$py_runner" -m pip --version; then
    echo "INFO: Installing pip for the selected Python."
    "$py_runner" -m ensurepip --upgrade || exit $?
fi

echo "INFO: Installing or updating EXOTIC."
"$py_runner" -m pip install --upgrade exotic || exit $?

echo "INFO: Launching EXOTIC user interface."
# Use the installed package rather than a checkout in the current directory.
exec "$py_runner" -P -c 'from exotic.exotic_gui import main; main()'
