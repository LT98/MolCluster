"""Keep the experiment log — see `mofsbu.explog.cli` for the commands.

    python scripts/expt.py check
    python scripts/expt.py new 2-04 --like 2-01

Runs this checkout's code against this checkout's data root, with no install needed.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.explog.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
