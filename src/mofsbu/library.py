"""A named molecule library — quick reference for things you use repeatedly.

Display and convenience only, exactly like fragment aliases: a library entry is a name
you gave a SMILES, and nothing in identity, search or the pipeline depends on one
existing.  Stored as JSON under `data/reference/` because that is the tracked part of
`data/`, so the library travels between machines while generated registries do not.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from mofsbu.config import REPO_ROOT

LIBRARY_PATH = REPO_ROOT / "data" / "reference" / "molecule_library.json"


@dataclass(frozen=True)
class Entry:
    name: str
    smiles: str
    multiplicity: int = 1
    note: str = ""
    added_at: str = ""


def load(path: Path | None = None) -> list[Entry]:
    p = Path(path or LIBRARY_PATH)
    if not p.exists():
        return []
    raw = json.loads(p.read_text(encoding="utf-8"))
    return [Entry(**e) for e in raw.get("molecules", [])]


def save(entries: list[Entry], path: Path | None = None) -> Path:
    p = Path(path or LIBRARY_PATH)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"molecules": [asdict(e) for e in entries]},
                            indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return p


def add(name: str, smiles: str, *, multiplicity: int = 1, note: str = "",
        path: Path | None = None) -> Entry:
    """Add or replace an entry.  The SMILES is parsed first: an unparseable one is refused."""
    from mofsbu.graph.from_mol import mol_from_smiles

    mol_from_smiles(smiles)                       # raises on invalid input
    entry = Entry(name=name.strip(), smiles=smiles.strip(), multiplicity=int(multiplicity),
                  note=note.strip(),
                  added_at=datetime.now(timezone.utc).isoformat(timespec="seconds"))
    entries = [e for e in load(path) if e.name != entry.name]
    entries.append(entry)
    entries.sort(key=lambda e: e.name.lower())
    save(entries, path)
    return entry


def remove(name: str, path: Path | None = None) -> bool:
    entries = load(path)
    kept = [e for e in entries if e.name != name]
    if len(kept) == len(entries):
        return False
    save(kept, path)
    return True
