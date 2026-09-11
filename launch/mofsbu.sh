#!/usr/bin/env bash
#
# Start the mofsbu viewer.  Double-click it, or run it from anywhere.
#
#   ./launch/mofsbu.sh                 CPU, finds the registry on its own
#   ./launch/mofsbu.sh --gpu           use the first CUDA device for ml_go runs
#   ./launch/mofsbu.sh --gpu --workers 8
#   ./launch/mofsbu.sh --db scratch.db --port 8001
#
# Everything here exists to replace three commands that each had a way to go wrong:
# `conda activate mofsbu` (the env is not always called that, and `activate` needs a
# shell that has been `conda init`-ed, which a double-click does not give you),
# `python scripts/viewer.py` (which python?), and remembering MOFSBU_DEVICE.
#
# The interpreter is used BY FULL PATH and never activated.  `conda activate` in a
# non-interactive shell fails in several different ways depending on how conda was
# installed, and none of the failures says so clearly.  `<env>/bin/python` needs no
# shell setup at all and behaves identically for our purposes.
set -euo pipefail

# ── where is the repo? ───────────────────────────────────────────────────────
# Resolve symlinks: this script is meant to be linked to from a desktop or a PATH
# directory, and $0 is then the link rather than the file.
SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do
  DIR=$(cd -P "$(dirname "$SOURCE")" && pwd)
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE
done
REPO=$(cd -P "$(dirname "$SOURCE")/.." && pwd)

# ── options ──────────────────────────────────────────────────────────────────
DEVICE=""; WORKERS=""; PORT="8000"; OPEN="--open"; PASS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --gpu)         DEVICE="cuda"; shift ;;
    --cpu)         DEVICE="cpu"; shift ;;
    --device)      DEVICE="$2"; shift 2 ;;
    --workers)     WORKERS="$2"; shift 2 ;;
    --port)        PORT="$2"; shift 2 ;;
    --no-browser)  OPEN=""; shift ;;
    -h|--help)
      sed -n '2,16p' "$SOURCE" | sed 's/^# \{0,1\}//'
      echo
      echo "Anything else is passed straight through to scripts/viewer.py:"
      echo "  --db FILE  --store DIR  --host ADDR  --list-db  --log-level LEVEL"
      exit 0 ;;
    *)             PASS+=("$1"); shift ;;
  esac
done

# ── which python? ────────────────────────────────────────────────────────────
# Probed, not assumed.  environment.yml says the env is called `mofsbu`, but an env is
# whatever it was actually created as — on at least one machine here it is not that —
# so the name is a hint and the IMPORTS are the test.  A wrong guess that starts and
# then fails at the first SMILES is worse than a slow, explicit search.
usable() { [ -x "$1" ] && "$1" -c 'import fastapi, uvicorn, rdkit' >/dev/null 2>&1; }

find_python() {
  # 1. An explicit choice always wins.
  if [ -n "${MOFSBU_PYTHON:-}" ]; then echo "$MOFSBU_PYTHON"; return; fi

  # 2. Named conda environments, under every root conda is normally installed to.
  local roots=() names=() root name
  [ -n "${CONDA_ROOT:-}" ] && roots+=("$CONDA_ROOT")
  roots+=("$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3"
          "/opt/miniforge3" "/opt/conda" "/usr/local/miniconda3")
  [ -n "${MOFSBU_ENV:-}" ] && names+=("$MOFSBU_ENV")
  names+=("mofsbu")
  # The name recorded in environment.yml, in case it is ever changed there.
  if [ -f "$REPO/environment.yml" ]; then
    name=$(sed -n 's/^name:[[:space:]]*//p' "$REPO/environment.yml" | head -1)
    [ -n "$name" ] && names+=("$name")
  fi
  for root in "${roots[@]}"; do
    for name in "${names[@]}"; do
      usable "$root/envs/$name/bin/python" && { echo "$root/envs/$name/bin/python"; return; }
    done
  done

  # 3. Any environment that can actually do the job.  This is the case that matters on a
  #    machine where the env was created under a different name.
  for root in "${roots[@]}"; do
    [ -d "$root/envs" ] || continue
    for candidate in "$root"/envs/*/bin/python; do
      usable "$candidate" && { echo "$candidate"; return; }
    done
  done

  # 4. Whatever is on PATH, including an already-active env.
  for candidate in python3 python; do
    command -v "$candidate" >/dev/null 2>&1 || continue
    usable "$(command -v "$candidate")" && { command -v "$candidate"; return; }
  done
}

PY=$(find_python || true)
if [ -z "${PY:-}" ]; then
  cat >&2 <<EOF

mofsbu: could not find a Python with the dependencies installed.

Looked for a conda environment containing fastapi, uvicorn and rdkit, and found none.
Create one once, from the repo root:

    conda env create -f "$REPO/environment.yml"

If the environment exists under a different name, say which:

    MOFSBU_ENV=the-name "$SOURCE"

or point straight at its interpreter:

    MOFSBU_PYTHON=/path/to/envs/NAME/bin/python "$SOURCE"

EOF
  # A double-clicked window closes the instant the script exits, taking the message with
  # it, so hold it open when there is nobody watching a terminal that will persist.
  [ -t 0 ] && [ -t 1 ] || { echo "Press Enter to close." >&2; read -r _ || true; }
  exit 1
fi

# ── go ───────────────────────────────────────────────────────────────────────
ARGS=(--port "$PORT")
[ -n "$OPEN" ] && ARGS+=("$OPEN")
[ -n "$DEVICE" ] && ARGS+=(--device "$DEVICE")
[ -n "$WORKERS" ] && ARGS+=(--workers "$WORKERS")
[ ${#PASS[@]} -gt 0 ] && ARGS+=("${PASS[@]}")

echo "mofsbu: $PY"
cd "$REPO"
# Unbuffered, because Python buffers stdout whenever it is not a terminal — and a
# launcher's stdout frequently is not (a file-manager window, a log file, a pipe).  The
# startup banner carries the URL, so it must not sit in a buffer until the server exits.
export PYTHONUNBUFFERED=1
exec "$PY" scripts/viewer.py "${ARGS[@]}"
