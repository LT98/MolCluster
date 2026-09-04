"""Run the reference cases in `data/reference/ligand_cases.tsv` and report.

    python scripts/check_cases.py                 # check every case
    python scripts/check_cases.py --name catechol # one of them, verbosely
    python scripts/check_cases.py --verbose       # show what was found, pass or fail

Exits non-zero if any case disagrees, so it can gate a commit.  Add cases by editing the
TSV; nothing here needs changing.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from mofsbu.geometry.embed import embed_molecule                           # noqa: E402
from mofsbu.graph.from_mol import mol_from_smiles                          # noqa: E402
from mofsbu.sites.model import chelate_pockets, perceive                   # noqa: E402
from mofsbu.sites.protomers import enumerate_protomers                     # noqa: E402

CASES = ROOT / "data" / "reference" / "ligand_cases.tsv"


def load_cases(path: Path = CASES) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n"))
    reader = csv.DictReader(rows, delimiter="\t")
    return [r for r in reader if r.get("smiles")]


def evaluate(case: dict) -> tuple[bool, str, str]:
    """Returns (ok, expected, found)."""
    mol = embed_molecule(mol_from_smiles(case["smiles"]), seed=7)
    sites = perceive(mol)
    pockets = chelate_pockets(mol, sites)
    protomers = enumerate_protomers(mol_from_smiles(case["smiles"]),
                                    multiplicity=int(case["multiplicity"]))

    problems: list[str] = []
    chelates = bool(pockets)
    if case["chelates"] in ("yes", "no"):
        want = case["chelates"] == "yes"
        if chelates is not want:
            problems.append(f"chelates={chelates} want={want}")
    if case["ring"] != "-" and pockets:
        rings = {p.ring_size for p in pockets}
        if int(case["ring"]) not in rings:
            problems.append(f"rings={sorted(rings)} want={case['ring']}")
    if case["donors"] != "-" and pockets:
        found = {"+".join(sorted(p.donor_types)) for p in pockets}
        if case["donors"] not in found:
            problems.append(f"donors={sorted(found)} want={case['donors']}")
    if case["protomers"] != "-":
        if len(protomers) != int(case["protomers"]):
            problems.append(f"protomers={len(protomers)} want={case['protomers']}")

    found_desc = (f"{len(pockets)} pocket(s) "
                  + ", ".join(sorted({p.descriptor for p in pockets}))
                  + f" | {len(protomers)} protomers")
    return not problems, "; ".join(problems), found_desc


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--name", help="run only this case")
    p.add_argument("--verbose", "-v", action="store_true")
    p.add_argument("--file", type=Path, default=CASES)
    a = p.parse_args(argv)

    cases = load_cases(a.file)
    if a.name:
        cases = [c for c in cases if c["name"] == a.name]
        if not cases:
            p.error(f"no case named {a.name!r} in {a.file}")

    failures = 0
    print(f"{a.file}  ({len(cases)} cases)\n")
    for case in cases:
        ok, problems, found = evaluate(case)
        mark = "ok  " if ok else "FAIL"
        print(f"  {mark} {case['name']:22s} {case['note']}")
        if not ok:
            failures += 1
            print(f"         -> {problems}")
            print(f"         found: {found}")
        elif a.verbose:
            print(f"         found: {found}")
    print(f"\n  {len(cases) - failures}/{len(cases)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
