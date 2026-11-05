"""A change to code must say what it did in the plan archive, in the same branch.

Compares the branch (committed and uncommitted) with where it left `origin/main`. If anything
under `src/` or `scripts/` changed, `docs/archive/PLAN_completed.md` must have changed too:
that is where the changelog revision for the work goes (AGENTS.md, "Finishing a change").

A change with nothing to record (a typo, a rename, a test-only fix) opts out with a commit
trailer that says why:

    Plan-update: none — <reason>

    python scripts/check_docs.py              # against origin/main
    python scripts/check_docs.py --base main
"""
from __future__ import annotations

import subprocess
import sys

CODE = ("src/", "scripts/")
CHANGELOG = "docs/archive/PLAN_completed.md"
OPT_OUT = "Plan-update: none"


def git(*args: str) -> str:
    return subprocess.run(("git", *args), capture_output=True, text=True,
                          check=True).stdout


def main(argv: list[str]) -> int:
    target = argv[argv.index("--base") + 1] if "--base" in argv else "origin/main"
    try:
        base = git("merge-base", "HEAD", target).strip()
    except subprocess.CalledProcessError:
        print(f"check_docs: no merge base with {target}; skipped")
        return 0
    changed = set(git("diff", "--name-only", base).split())
    changed |= set(git("ls-files", "--others", "--exclude-standard").split())
    code = sorted(p for p in changed if p.startswith(CODE))
    if not code or CHANGELOG in changed:
        return 0
    if OPT_OUT in git("log", "--format=%B", f"{base}..HEAD"):
        return 0
    print(f"check_docs: this branch changes {len(code)} code file(s) "
          f"({', '.join(code[:4])}{', ...' if len(code) > 4 else ''})\n"
          f"but not {CHANGELOG}. Add a changelog revision for the work, and update\n"
          f"PLAN_implementation.md / BUGS / DESIGN / the WORKPLAN it advances (AGENTS.md,\n"
          f"'Finishing a change'). Nothing to record? Add a commit trailer:\n"
          f"    {OPT_OUT} — <reason>")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
