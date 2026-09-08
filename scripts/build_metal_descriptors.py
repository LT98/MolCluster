#!/usr/bin/env python3
"""Generate `data/reference/metal_descriptors.tsv` (M1, decision C8).

The split this implements: **ionic radii come from a package, everything else is
curated.**  Shannon radii are genuinely standard, genuinely tabulated, and genuinely
tedious to retype per charge / coordination number / spin state — exactly what a
reference package is for.  HSAB class, preferred coordination numbers and water-exchange
lability are NOT in any package, are chemical judgement, and are typed out below with
their reasoning visible.

The generated TSV is what ships; `mendeleev` is a dev dependency, not a runtime one, so a
laptop with no scientific stack still gets the table by cloning the repo.  Re-running this
with a different `mendeleev` bumps `source_version` in the output, which turns "the
numbers moved" into a visible diff rather than a silent one (ground rule 6).

    pip install mendeleev
    python scripts/build_metal_descriptors.py
"""
from __future__ import annotations

import argparse
import importlib.metadata
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

OUT = ROOT / "data" / "reference" / "metal_descriptors.tsv"

# symbol, charge -> (hsab, preferred_cn, exchange_lability, spin for the radius lookup)
#
# hsab follows Pearson's original classification.  exchange_lability is the water-exchange
# rate class: fast ~1e8-1e9 /s (Cu2+ Jahn-Teller, Zn2+, Ca2+), moderate ~1e4-1e7,
# slow <1e2 (Cr3+ at ~1e-6 is the textbook inert ion).  It is a C6 barrier-proxy input:
# a slow ion pays a desolvation penalty a thermodynamic dG never sees.
CURATED: dict[tuple[str, int], tuple[str, tuple[int, ...], str, str]] = {
    ("Mg", 2): ("hard",       (6,),      "moderate", ""),
    ("Ca", 2): ("hard",       (6, 7, 8), "fast",     ""),
    ("Al", 3): ("hard",       (4, 6),    "slow",     ""),
    ("Sc", 3): ("hard",       (6,),      "fast",     ""),
    ("Ti", 4): ("hard",       (6,),      "slow",     ""),
    ("V",  3): ("hard",       (6,),      "moderate", ""),
    ("V",  4): ("hard",       (5, 6),    "moderate", ""),
    ("Cr", 3): ("hard",       (6,),      "slow",     ""),
    ("Mn", 2): ("borderline", (6,),      "fast",     "HS"),
    ("Fe", 2): ("borderline", (6,),      "moderate", "HS"),
    ("Fe", 3): ("hard",       (6,),      "moderate", "HS"),
    ("Co", 2): ("borderline", (4, 6),    "moderate", "HS"),
    ("Co", 3): ("hard",       (6,),      "slow",     "LS"),
    ("Ni", 2): ("borderline", (4, 6),    "moderate", ""),
    ("Cu", 1): ("soft",       (2, 4),    "fast",     ""),
    ("Cu", 2): ("borderline", (4, 5, 6), "fast",     ""),
    ("Zn", 2): ("borderline", (4, 6),    "fast",     ""),
    ("Ga", 3): ("hard",       (4, 6),    "moderate", ""),
    ("Y",  3): ("hard",       (8,),      "fast",     ""),
    ("Zr", 4): ("hard",       (7, 8),    "slow",     ""),
    ("Ru", 3): ("borderline", (6,),      "slow",     ""),
    ("Rh", 3): ("borderline", (6,),      "slow",     ""),
    ("Pd", 2): ("soft",       (4,),      "slow",     ""),
    ("Ag", 1): ("soft",       (2, 4),    "fast",     ""),
    ("Cd", 2): ("soft",       (6,),      "fast",     ""),
    ("In", 3): ("hard",       (6,),      "moderate", ""),
    ("La", 3): ("hard",       (8, 9),    "fast",     ""),
    ("Ce", 3): ("hard",       (8, 9),    "fast",     ""),
    ("Ce", 4): ("hard",       (8,),      "moderate", ""),
    ("Hf", 4): ("hard",       (7, 8),    "slow",     ""),
    ("Pt", 2): ("soft",       (4,),      "slow",     ""),
    ("Au", 1): ("soft",       (2,),      "fast",     ""),
    ("Hg", 2): ("soft",       (2, 4),    "fast",     ""),
    ("Pb", 2): ("borderline", (6,),      "fast",     ""),
    ("Bi", 3): ("borderline", (6, 8),    "fast",     ""),
}


def shannon_radius(symbol: str, charge: int, cn: int, spin: str) -> float | None:
    """Shannon effective ionic radius, in angstrom, at this CN and spin.

    Falls back across coordination numbers rather than inventing a number: an ion with no
    tabulated radius at its preferred CN returns the nearest tabulated one, and one with
    nothing at all returns None so the row is visibly incomplete.
    """
    from mendeleev import element

    entries = [r for r in element(symbol).ionic_radii if r.charge == charge]
    if not entries:
        return None
    if spin:
        entries = [r for r in entries if (r.spin or "") == spin] or entries
    exact = [r for r in entries if _cn_of(r) == cn]
    pool = exact or entries
    # Shannon flags the values he considered well determined.  Prefer those, then the
    # coordination number actually asked for, then the nearest one -- and never average
    # across coordination numbers, which would invent a radius no table contains.
    reliable = [r for r in pool if getattr(r, "most_reliable", False)]
    pool = reliable or pool
    chosen = pool[0] if exact else min(pool, key=lambda r: abs((_cn_of(r) or 99) - cn))
    return round(chosen.ionic_radius / 100.0, 3)      # pm -> angstrom


_ROMAN = {"II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6, "VII": 7, "VIII": 8,
          "IX": 9, "X": 10, "XI": 11, "XII": 12}


def _cn_of(entry) -> int:
    raw = str(entry.coordination or "").strip().rstrip("SPY")
    return _ROMAN.get(raw, 0)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=OUT)
    args = p.parse_args(argv)

    try:
        version = importlib.metadata.version("mendeleev")
    except importlib.metadata.PackageNotFoundError:
        raise SystemExit(
            "this generator needs `pip install mendeleev` (a dev dependency; the "
            "generated TSV is what ships)") from None

    from mofsbu.energy.backends import d_electrons

    lines = [
        "# metal_descriptors — GENERATED, do not hand-edit.",
        f"#   python scripts/build_metal_descriptors.py     (mendeleev {version})",
        "#",
        "# ionic_radius: Shannon effective ionic radius in angstrom, at the FIRST",
        "#   preferred_cn and, where it matters, the high-spin state -- from mendeleev.",
        "# hsab / preferred_cn / exchange_lability: curated in the generator, because no",
        "#   package carries them.  See the CURATED table there for the reasoning.",
        "# d_electrons: derived from group minus charge; empty for main-group ions.",
        "symbol\tcharge\tionic_radius\thsab\tpreferred_cn\td_electrons\texchange_lability\tsource\tsource_version",
    ]
    missing: list[str] = []
    for (symbol, charge), (hsab, cns, lability, spin) in sorted(CURATED.items()):
        radius = shannon_radius(symbol, charge, cns[0], spin)
        if radius is None:
            missing.append(f"{symbol}{charge:+d}")
        try:
            d = d_electrons(symbol, charge)
        except ValueError:
            d = None                                   # main group: no d-count to give
        lines.append("\t".join([
            symbol, str(charge), "" if radius is None else f"{radius:.3f}", hsab,
            ",".join(str(c) for c in cns), "" if d is None else str(d), lability,
            "mendeleev+curated", version,
        ]))

    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {args.out}  ({len(CURATED)} ions, mendeleev {version})")
    if missing:
        print("no tabulated Shannon radius for: " + ", ".join(missing)
              + "  — those rows carry an empty radius rather than a guess")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
