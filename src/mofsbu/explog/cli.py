"""`expt` — keep the experiment log.

    python scripts/expt.py init                         # folders, reagents.yaml, phases.yaml
    python scripts/expt.py new 2-04 --like 2-01         # copy the recipe, blank the outcome
    python scripts/expt.py new 4-01 --method solution   # blank template for a method
    python scripts/expt.py fmt [ID ...]                 # rewrite in canonical, commented form
    python scripts/expt.py check [ID ...]               # what is still to record; problems
    python scripts/expt.py show 2-01                    # one entry as a readable card
    python scripts/expt.py diff 2-01 2-02               # only the fields that differ
    python scripts/expt.py index                        # regenerate INDEX.md

The log lives in `<data root>/experiments` (`--root` to point elsewhere).  Every command
refreshes `.schema/` so the editor's dropdowns follow reagents.yaml and phases.yaml.
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

from ..config import data_root
from .emit import render_entry, render_phases, render_reagents
from .model import Log, check_entry, check_phases, check_reagents, load_log
from .render import diff, index, show
from .schema import write_schemas
from .spec import METHODS, OUTCOME_SECTIONS


class UsageError(Exception):
    pass


def _ids(log: Log, ids: list[str]) -> list[str]:
    unknown = [i for i in ids if i not in log.entries]
    if unknown:
        raise UsageError(f"no entry {', '.join(unknown)} under {log.root / 'entries'}")
    return ids or list(log.entries)


def _refuse_load_errors(log: Log) -> None:
    if log.load_errors:
        raise UsageError("cannot read: " + "; ".join(f"{n}: {e}"
                                                      for n, e in log.load_errors.items()))


def cmd_init(log: Log, args) -> int:
    for sub in ("entries", "raw"):
        (log.root / sub).mkdir(parents=True, exist_ok=True)
    for name, render in (("reagents", render_reagents), ("phases", render_phases)):
        path = log.root / f"{name}.yaml"
        if not path.exists():
            path.write_text(render({}), encoding="utf-8")
            print(f"wrote {path}")
    print(f"log ready at {log.root}")
    return 0


def cmd_new(log: Log, args) -> int:
    path = log.entry_path(args.id)
    if path.exists():
        raise UsageError(f"{path} already exists")
    if args.like:
        if args.like not in log.entries:
            raise UsageError(f"no entry {args.like} to copy")
        values = copy.deepcopy(log.entries[args.like])
        for key in ("title", "date", "notebook_ref", "notes", *OUTCOME_SECTIONS):
            values[key] = None
        values.update(id=args.id, derived_from=args.like, status="planned")
    else:
        if not args.method:
            raise UsageError("give --method or --like")
        values = {"id": args.id, "status": "planned", "method": args.method}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_entry(values), encoding="utf-8")
    print(f"wrote {path}")
    return 0


def cmd_fmt(log: Log, args) -> int:
    _refuse_load_errors(log)
    targets = [(log.files[i], render_entry(log.entries[i])) for i in _ids(log, args.ids)]
    if not args.ids:
        for name, table, render in (("reagents", log.reagents, render_reagents),
                                    ("phases", log.phases, render_phases)):
            if (log.root / f"{name}.yaml").exists():
                targets.append((log.root / f"{name}.yaml", render(table)))
    for path, text in targets:
        if path.read_text(encoding="utf-8") != text:
            path.write_text(text, encoding="utf-8")
            print(f"formatted {path.relative_to(log.root)}")
    return 0


def cmd_check(log: Log, args) -> int:
    bad = bool(log.load_errors)
    for name, err in log.load_errors.items():
        print(f"{name}: cannot read: {err}")
    ids = _ids(log, args.ids)
    groups = [(i, check_entry(log, i)) for i in ids]
    if not args.ids:
        groups += [(f"reagents.yaml {c}", f) for c, f in check_reagents(log).items()]
        groups += [(f"phases.yaml {c}", f) for c, f in check_phases(log).items()]
    for name, found in groups:
        if not (found.problems or found.pending and not args.problems_only):
            continue
        print(f"\n{name}: {len(found.pending)} to record, {len(found.problems)} problems")
        for p in found.problems:
            print(f"  ! {p}")
        if not args.problems_only:
            for p in found.pending:
                print(f"  - {p}")
        bad = bad or bool(found.problems)
    return 1 if bad else 0


def cmd_show(log: Log, args) -> int:
    print(show(log, _ids(log, [args.id])[0]), end="")
    return 0


def cmd_diff(log: Log, args) -> int:
    _ids(log, [args.a, args.b])
    print(diff(log, args.a, args.b), end="")
    return 0


def cmd_index(log: Log, args) -> int:
    path = log.root / "INDEX.md"
    path.write_text(index(log), encoding="utf-8")
    print(f"wrote {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="expt", description=__doc__.split("\n")[0])
    parser.add_argument("--root", type=Path, help="log folder (default <data root>/experiments)")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("init", help="create the folder layout")
    p = sub.add_parser("new", help="start an entry")
    p.add_argument("id")
    p.add_argument("--like", help="copy the recipe of this entry")
    p.add_argument("--method", choices=METHODS)
    p = sub.add_parser("fmt", help="rewrite entries in canonical form")
    p.add_argument("ids", nargs="*")
    p = sub.add_parser("check", help="list what is still to record, and problems")
    p.add_argument("ids", nargs="*")
    p.add_argument("--problems-only", action="store_true")
    p = sub.add_parser("show", help="print one entry")
    p.add_argument("id")
    p = sub.add_parser("diff", help="print fields that differ between two entries")
    p.add_argument("a")
    p.add_argument("b")
    sub.add_parser("index", help="regenerate INDEX.md")
    args = parser.parse_args(argv)

    root = args.root or data_root() / "experiments"
    log = load_log(root)
    if args.command != "init" and not root.exists():
        print(f"no log at {root}; run `expt init` first", file=sys.stderr)
        return 2
    handlers = {"init": cmd_init, "new": cmd_new, "fmt": cmd_fmt, "check": cmd_check,
                "show": cmd_show, "diff": cmd_diff, "index": cmd_index}
    try:
        status = handlers[args.command](log, args)
    except UsageError as exc:
        print(f"expt: {exc}", file=sys.stderr)
        return 2
    write_schemas(load_log(root))
    return status
