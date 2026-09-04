"""Content-addressed blob store (D2: address is a content hash, meaning lives in DB rows).

Bytes go in, a sha256 digest comes out, and the digest is the only handle.  Filenames
carry no meaning, so the length limit that shaped the legacy pipeline cannot come back.
Writes are atomic and idempotent: storing identical bytes twice is one file.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

SHARD = 2          # store/ab/cdef… — keeps directory sizes sane


class BlobStore:
    __slots__ = ("root",)

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def digest(data: bytes) -> str:
        return hashlib.sha256(data).hexdigest()

    def path(self, digest: str) -> Path:
        if len(digest) != 64 or not all(c in "0123456789abcdef" for c in digest):
            raise ValueError(f"not a sha256 hex digest: {digest!r}")
        return self.root / digest[:SHARD] / f"{digest[SHARD:]}.blob"

    def put(self, data: bytes) -> str:
        d = self.digest(data)
        p = self.path(d)
        if p.exists():
            return d                                  # idempotent
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, p)                            # atomic
        return d

    def put_text(self, text: str) -> str:
        return self.put(text.encode("utf-8"))

    def get(self, digest: str) -> bytes:
        p = self.path(digest)
        if not p.exists():
            raise KeyError(f"blob {digest} not in {self.root}")
        return p.read_bytes()

    def get_text(self, digest: str) -> str:
        return self.get(digest).decode("utf-8")

    def has(self, digest: str) -> bool:
        return self.path(digest).exists()

    def __len__(self) -> int:
        return sum(1 for _ in self.root.glob(f"{'?' * SHARD}/*.blob"))
