#!/usr/bin/env bash

# Use the active environment and run this checkout directly.
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
cd -- "$(dirname -- "$0")" || exit $?
"$py_runner" -m exotic.api.native_build || exit $?
echo "INFO: Launching EXOTIC user interface from this checkout."
exec "$py_runner" -c 'from exotic.exotic_gui import main; main()'
