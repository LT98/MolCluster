#!/usr/bin/env python3
"""M7 exit gate — does the new energy stack reproduce the archived work?

`docs/PLAN_implementation.md` §M7 sets the gate: re-run the Ni/BTC and Fe/BTC cases
through the new stack and reproduce the rankings in `docs/reports/*.md` within noise.
"If the new stack can't reproduce the old results, one of them is wrong and you want to
know NOW, not in M8."

Run from the repo root, on a machine with tblite:

    python scripts/regress_m7.py --refs           # stage 1, minutes
    python scripts/regress_m7.py --refs --json out.json

Stage 1 is the cheap, exact half.  `legacy/fe_btc_refs.json` records the four reference
energies the archived Fe(III) run subtracted — the sextet Fe(3+) ion and the BTC, EDTA
and EDDA anions.  Same method, same species: recomputing them is a direct check that the
new `XTBBackend` is the same calculator that produced the archived numbers.  It is not a
check of the ranking, and it does not pretend to be.

The comparison is deliberately not bit-exact.  The archived ligand geometries came from
`ebu_core.LigandBuilder`; these come from `mofsbu.geometry.embed`.  Two ETKDG embeddings
of EDTA are two different conformers, and a light LBFGS does not erase that, so a few
tenths of an eV of disagreement is the embedding and not the energy model.  A whole eV is
not, and that is what the tolerance is set to catch.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from mofsbu._types import EnergyBackendUnavailable                  # noqa: E402
from mofsbu.energy import XTBBackend, high_spin_multiplicity        # noqa: E402
from mofsbu.geometry.embed import embed_molecule                    # noqa: E402
from mofsbu.graph.from_mol import mol_from_smiles                   # noqa: E402

# The archived run's ligands, as deprotonated as the reference file assumed: the anion,
# not the acid.  SMILES carry the charge explicitly rather than relying on a pH rule.
LIGANDS = {
    "BTC":  "[O-]C(=O)c1cc(C(=O)[O-])cc(C(=O)[O-])c1",
    "EDTA": "C(CN(CC(=O)[O-])CC(=O)[O-])N(CC(=O)[O-])CC(=O)[O-]",
    "EDDA": "[O-]C(=O)CNCCNCC(=O)[O-]",
}
ARCHIVED = ROOT / "legacy" / "fe_btc_refs.json"

# The Fe(III) run fixed a sextet for every Fe species (docs/reports/fe_btc_report.md).
# That convention is now `high_spin_multiplicity('Fe', 3)` rather than an env var.
METAL, METAL_CHARGE = "Fe", 3


def _relax(backend, symbols, positions, *, charge, multiplicity, fmax, steps):
    return backend.relax(symbols, positions, charge=charge, multiplicity=multiplicity,
                         fmax=fmax, steps=steps)


def stage_refs(*, tol: float, fmax: float, steps: int) -> dict:
    if not ARCHIVED.exists():
        raise SystemExit(f"no archived reference file at {ARCHIVED}")
    archived = json.loads(ARCHIVED.read_text())

    backend = XTBBackend()
    if not backend.available():
        raise EnergyBackendUnavailable(
            f"stage 1 needs GFN2-xTB: {backend.install_hint()}")

    mult = high_spin_multiplicity(METAL, METAL_CHARGE)
    print(f"# GFN2-xTB via {backend.code}-{backend.code_version()}")
    print(f"# {METAL}({METAL_CHARGE}+) as multiplicity {mult} "
          f"(high spin d{10 - METAL_CHARGE - 2}, the archived sextet convention)\n")
    print(f"{'species':10s} {'archived':>12s} {'recomputed':>12s} {'delta':>9s}  verdict")

    rows: list[dict] = []
    ion = backend.single_point([METAL], [[0.0, 0.0, 0.0]],
                               charge=METAL_CHARGE, multiplicity=mult)
    rows.append({"species": METAL, "archived": archived["M"], "recomputed": ion.energy})

    for name, smiles in LIGANDS.items():
        mol = embed_molecule(mol_from_smiles(smiles), seed=7)
        conf = mol.GetConformer()
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        positions = conf.GetPositions()
        charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        result = _relax(backend, symbols, positions, charge=charge, multiplicity=1,
                        fmax=fmax, steps=steps)
        rows.append({"species": name, "archived": archived["lig"][name],
                     "recomputed": result.energy, "charge": charge,
                     "converged": result.converged, "steps": result.n_steps})

    worst = 0.0
    for row in rows:
        delta = row["recomputed"] - row["archived"]
        row["delta"] = delta
        row["within_tolerance"] = abs(delta) <= tol
        worst = max(worst, abs(delta))
        mark = "ok" if row["within_tolerance"] else "OUT OF TOLERANCE"
        print(f"{row['species']:10s} {row['archived']:12.3f} {row['recomputed']:12.3f} "
              f"{delta:+9.3f}  {mark}")

    passed = all(r["within_tolerance"] for r in rows)
    print(f"\nworst |delta| = {worst:.3f} eV against a tolerance of {tol:.3f} eV")
    if passed:
        print("STAGE 1 PASS — the new backend reproduces the archived reference energies.")
    else:
        print("STAGE 1 FAIL — the same species at the same level of theory disagree by "
              "more than a conformer can explain. Before trusting either set of numbers, "
              "check the charge and multiplicity on the row that failed: those are the "
              "two things the archived run passed through environment variables.")
    return {"stage": "refs", "tolerance": tol, "passed": passed,
            "multiplicity": mult, "rows": rows,
            "code_version": backend.code_version()}


def stage_rankings() -> dict:
    """Stage 2: reproduce the Ni/BTC and Fe/BTC landscape rankings.  NOT BUILT.

    Ground rule 8 — the signature is settled, the body is scheduled work, and it does not
    return a partial answer dressed as a whole one.  What it needs that does not exist:

    * the archived candidate geometries.  `docs/reports/ni_btc_report.md` lists
      `ni_btc_outputs/` and `ni_btc_xtb_results.csv` as its files; neither is in the repo,
      and `legacy/ni_btc_results.csv` has an EMPTY `formation_eV` column.  The ranking to
      regress against currently exists only as the table in the report.
    * a re-enumeration through the new stack to replace them, which is a `BuildSpec` run
      per scenario — cheap to write and hours of xTB to execute, so it belongs behind its
      own flag and its own resumable task table, not inside this function.
    * the reference scheme to score them with.  `energy.reference` will refuse the
      archived E(EBU) - E(M) - SUM E(anion) equation as non-isodesmic, which is correct
      and which means stage 2 is not a re-run: it is a re-derivation under a different
      scheme, and the honest comparison is of ORDERINGS, not of values.
    """
    from mofsbu.assembly.join import NotBuiltYet

    raise NotBuiltYet(
        "regress_m7 stage 2 (landscape rankings) — needs the archived candidate "
        "geometries, which are not in the repo, or a re-enumeration to replace them. "
        "Run stage 1 (--refs) first; it is the half that can be checked exactly.")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refs", action="store_true",
                   help="stage 1: reproduce the archived Fe(III) reference energies")
    p.add_argument("--rankings", action="store_true",
                   help="stage 2: reproduce the landscape rankings (not built — says why)")
    p.add_argument("--tol", type=float, default=1.0,
                   help="eV tolerance per species (default 1.0; a different ETKDG "
                        "conformer is worth a few tenths, a wrong charge is worth tens)")
    p.add_argument("--fmax", type=float, default=0.20, help="LBFGS force threshold, eV/A")
    p.add_argument("--steps", type=int, default=35,
                   help="LBFGS step cap (35 = the archived run's cap)")
    p.add_argument("--json", type=Path, default=None, help="also write the result here")
    args = p.parse_args(argv)

    if not (args.refs or args.rankings):
        p.error("choose --refs and/or --rankings")

    out: dict = {}
    if args.refs:
        out["refs"] = stage_refs(tol=args.tol, fmax=args.fmax, steps=args.steps)
    if args.rankings:
        stage_rankings()

    if args.json:
        args.json.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
        print(f"\nwrote {args.json}")
    return 0 if all(v.get("passed") for v in out.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
