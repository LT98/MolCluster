"""Each subpackage OVERVIEW.md's generated API section matches the code."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import gen_overview  # noqa: E402


def test_generated_overviews_are_current():
    assert gen_overview.main(["--check"]) == 0, \
        "an OVERVIEW.md is stale: run `python scripts/gen_overview.py` and commit the result"
