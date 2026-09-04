"""Read-only web viewer over the structure registry (M3.5).

The viewer is *output only*: it opens the database with a read-only SQLite URI and
never writes.  That is the whole risk story — it cannot corrupt the registry, so it
needs no job state and stays disposable (docs/PLAN_implementation.md §2.5).

    python -m mofsbu.ui --db data/demo_registry.db --store data/store
"""
from __future__ import annotations

from .app import create_app

__all__ = ["create_app"]
