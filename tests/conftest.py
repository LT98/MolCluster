from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:          # works without an editable install
    sys.path.insert(0, str(ROOT / "src"))

import fixtures as fx  # noqa: E402


@pytest.fixture(scope="session")
def all_fixtures():
    return {name: build() for name, build in fx.ALL.items()}


@pytest.fixture(scope="session")
def golden_path() -> Path:
    return Path(__file__).parent / "golden_hashes.json"
