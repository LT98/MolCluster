"""Start the registry viewer from a checkout, with or without an editable install.

    python scripts/viewer.py                       # defaults to the demo registry
    python scripts/viewer.py --db data/registry.db --port 8000

Same launcher pattern as the other scripts here: it puts `src` on the path itself, so
it works from a fresh clone before `pip install -e .` has been run.  After the install,
`python -m mofsbu.ui` and the `mofsbu-ui` command do the same thing.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.ui.__main__ import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
