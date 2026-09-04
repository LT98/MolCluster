"""The reference cases run as part of the suite, so adding one guards it permanently."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from check_cases import evaluate, load_cases  # noqa: E402

CASES = load_cases()


def test_the_reference_file_is_not_empty():
    assert len(CASES) >= 10


@pytest.mark.parametrize("case", CASES, ids=[c["name"] for c in CASES])
def test_reference_case(case):
    ok, problems, found = evaluate(case)
    assert ok, f"{case['name']}: {problems}\n  found: {found}\n  ({case['note']})"
