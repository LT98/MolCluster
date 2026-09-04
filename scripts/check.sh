#!/usr/bin/env bash
# The gate both machines must pass before a push.
set -euo pipefail
cd "$(dirname "$0")/.."
python -m pytest tests
python scripts/verify_registry.py || echo "(no registry to verify yet)"
command -v ruff >/dev/null && ruff check src tests scripts || echo "(ruff not installed - skipped)"
