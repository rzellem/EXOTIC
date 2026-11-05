#!/usr/bin/env bash
# ==============================================================================
#  EXOTIC GPU — one-shot CachyOS / Arch setup script
# ==============================================================================
#  Run from the repo root (the folder that contains EXOTIC/):
#    cd "/home/manuelr/Documentos/GitProjects/EXOTIC GPU"
#    bash EXOTIC/setup_cachyos_venv.sh
#
#  CRITICAL RULE for CachyOS baked-in NVIDIA kernels:
#    If pacman ever prompts about package conflicts between nvidia-open-dkms
#    and linux-cachyos-*-nvidia-open:
#    → press N (the default). NEVER install nvidia-dkms on these kernels.
# ==============================================================================
set -u
set -o pipefail

RED=$'\033[0;31m'; GRN=$'\033[0;32m'; YLW=$'\033[1;33m'; RST=$'\033[0m'
p_ok()  { printf "${GRN}[OK]${RST}  %s\n" "$*"; }
p_warn(){ printf "${YLW}[!!]${RST} %s\n" "$*"; }
p_err() { printf "${RED}[FAIL]${RST} %s\n" "$*" >&2; }

cd "$(dirname "$0")/.." || exit 1  # repo root (parent of EXOTIC/)
p_ok "Working directory: $(pwd)"

# ---------------------------------------------------------------- 1) Python (use the preinstalled system python3)
PY=python3
if ! command -v "$PY" >/dev/null 2>&1; then
  p_err "$PY not found. Install 'python' via pacman first."
  exit 2
fi
p_ok "Python binary: $(command -v "$PY")"
p_ok "Version: $("$PY" --version 2>&1)"

# Make sure python-pip + python-virtualenv system packages are present
# (so "python -m venv" and "pip" are guaranteed to work; the "python"
# package itself often already ships venv, so --needed = no-op if OK).
sudo pacman -Syu --needed --noconfirm python-pip python-virtualenv || true

# ---------------------------------------------------------------- 2) Venv
VENV_DIR="$PWD/.venv-exotic"
if [ ! -x "$VENV_DIR/bin/python" ]; then
  p_warn "Creating venv at $VENV_DIR"
  "$PY" -m venv "$VENV_DIR" || { p_err "venv creation failed"; exit 3; }
fi
p_ok "Venv exists: $VENV_DIR"

# Activate (strict — we never touch system python)
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate" || { p_err "could not source venv"; exit 4; }
p_ok "Venv activated — pip is: $(command -v pip)"
command -v pip | grep -q "venv-exotic/bin/pip" || {
  p_err "pip is NOT inside the venv — aborting before we touch system Python"
  exit 5
}
p_ok "Python in venv: $(python --version 2>&1)"

# ---------------------------------------------------------------- 3) Bootstrap pip / wheel
p_ok "Upgrading pip / setuptools / wheel (inside venv)"
pip install -U pip setuptools wheel || exit 6

# ---------------------------------------------------------------- 4) Core requirements
REQ_FILE="$PWD/EXOTIC/requirements.txt"
if [ -f "$REQ_FILE" ]; then
  p_ok "Installing EXOTIC requirements from $REQ_FILE"
  pip install -r "$REQ_FILE" || {
    p_err "'pip install -r requirements.txt' failed. Stderr above."
    exit 7
  }
else
  p_err "requirements.txt not found at $REQ_FILE — are you running from the right folder?"
  exit 8
fi

# ---------------------------------------------------------------- 5) PyTorch CUDA
p_ok "Installing PyTorch CUDA 12.1 (recommended backend for alignment + photometry)"
pip install torch --index-url https://download.pytorch.org/whl/cu121 || {
  p_err "PyTorch CUDA wheel failed to install"
  p_warn "Fallback: pip install torch --index-url https://download.pytorch.org/whl/cpu"
  exit 9
}

# ---------------------------------------------------------------- 6) NVIDIA driver status
echo
if command -v nvidia-smi >/dev/null 2>&1; then
  p_ok "nvidia-smi found — first 12 lines:"
  nvidia-smi | head -n 12 || true
else
  p_warn "nvidia-smi is missing. Is the NVIDIA driver package installed (linux-cachyos-*-nvidia-open)?"
fi

# ---------------------------------------------------------------- 7) Validator
echo
p_ok "Running EXOTIC/check_gpu_install.py — this is THE source of truth:"
python EXOTIC/check_gpu_install.py
STATUS=$?
echo
if [ "$STATUS" -eq 0 ]; then
  p_ok "Validator passed. Now run EXOTIC normally."
  echo "  Reminder — enable GPU:  export EXOTIC_ENABLE_GPU=1  (default is 1)"
  echo "             disable:    export EXOTIC_ENABLE_GPU=0"
else
  p_warn "Validator exit code $STATUS. Fix the FAIL / WARN sections above,"
  p_warn "or paste the full output back to the EXOTIC GPU assistant."
  exit 10
fi
