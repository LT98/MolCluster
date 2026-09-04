"""Test fixtures — the example structures now live in the package.

Kept as a shim so tests keep importing `fixtures`, while scripts and the demo
registry draw the same graphs from `mofsbu.examples`.
"""
from __future__ import annotations

from mofsbu.examples import *  # noqa: F401,F403
from mofsbu.examples import ALL  # noqa: F401
