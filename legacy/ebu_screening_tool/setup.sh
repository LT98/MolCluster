#!/usr/bin/env bash
# One-shot setup for a fresh Linux machine.
set -euo pipefail
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
echo "[1/3] creating virtual environment (.venv) ..."
"$PY" -m venv .venv
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --upgrade pip wheel >/dev/null
echo "[2/3] installing requirements (this pulls rdkit, tblite, ase, jupyter ...) ..."
pip install -r requirements.txt
echo "[3/3] running self-test ..."
python selftest.py
cat << 'MSG'

============================================================
 Setup complete.  To start working:

   source .venv/bin/activate
   jupyter lab ebu_tool.ipynb

 Or from Python / scripts (run from this folder):
   python run_metal.py           # enumerate a metal's EBUs
   python xtb_energy_metal.py     # rank formation energies
============================================================
MSG
