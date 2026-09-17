from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:          # works without an editable install
    sys.path.insert(0, str(ROOT / "src"))

import fixtures as fx  # noqa: E402


@pytest.fixture(autouse=True)
def _declarations_do_not_outlive_their_test():
    """Every `MOFSBU_*` declaration is undone at the end of the test that made it.

    Not tidiness — isolation.  These variables are how the device and the worker count
    reach a run (`config`), and several tests set them through the code under test rather
    than through monkeypatch: `POST /api/compute` writes `MOFSBU_WORKERS` into this
    process on purpose, because that is how the page reaches the work.  One such test
    left the whole session declared a workstation, so a LATER test that called `run()`
    spawned workers it never asked for — and a spawned worker does not inherit
    monkeypatched module state, so tests that serve a fake energy backend to an
    in-process worker silently got the real one and hung on it.

    That failure depended on file order, which is the worst property a test failure can
    have.  It is fixed here, once, rather than in each test that declares something.
    """
    import os

    before = {k: v for k, v in os.environ.items() if k.startswith("MOFSBU_")}
    yield
    for key in [k for k in os.environ if k.startswith("MOFSBU_")]:
        if key not in before:
            del os.environ[key]
    os.environ.update(before)


@pytest.fixture(scope="session")
def all_fixtures():
    return {name: build() for name, build in fx.ALL.items()}


@pytest.fixture(scope="session")
def golden_path() -> Path:
    return Path(__file__).parent / "golden_hashes.json"
