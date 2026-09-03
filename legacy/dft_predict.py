"""
DFT Energy Prediction from .xyz files using MACE-MP-0
======================================================
MACE-MP-0 is a universal equivariant GNN (Batatia et al. 2023) pretrained on
the Materials Project database (~150k structures). It handles transition metals
(Zn, Fe, etc.) and organic ligands — making it well-suited for MOF SBUs.

Install:
    pip install mace-torch ase

Usage:
    python dft_predict.py molecule.xyz
    python dft_predict.py MOF_SBU/ebu_outputs/*.xyz --forces
    python dft_predict.py molecule.xyz --model small  # small | medium | large
"""

import argparse
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Core prediction
# ---------------------------------------------------------------------------

def load_model(model_size: str = "small"):
    """Load MACE-MP-0 calculator. Downloads weights on first call (~50–200 MB)."""
    try:
        from mace.calculators import mace_mp
    except ImportError:
        sys.exit(
            "mace-torch not found.\n"
            "Install with:  pip install mace-torch ase"
        )

    # dispersion=True adds D3 correction (requires torch-dftd); fall back silently
    try:
        calc = mace_mp(model=model_size, dispersion=True, default_dtype="float64")
    except RuntimeError:
        calc = mace_mp(model=model_size, dispersion=False, default_dtype="float64")
    return calc


def read_xyz(path: str | Path):
    """Parse .xyz file into an ASE Atoms object."""
    try:
        from ase.io import read
    except ImportError:
        sys.exit("ase not found.\nInstall with:  pip install ase")

    atoms = read(str(path), format="xyz")
    return atoms


def predict(xyz_path: str | Path, calc, compute_forces: bool = False) -> dict:
    """
    Run MACE-MP-0 inference on a single .xyz file.

    Returns
    -------
    dict with keys:
        file        : str
        n_atoms     : int
        energy_eV   : float   — total potential energy in eV
        energy_Ha   : float   — same in Hartree
        energy_per_atom_eV : float
        forces_eV_Ang : np.ndarray | None  — shape (N, 3), if requested
        elements    : list[str]
    """
    EV_TO_HARTREE = 0.036749405469679

    atoms = read_xyz(xyz_path)
    atoms.calc = calc

    energy_eV = float(atoms.get_potential_energy())

    forces = None
    if compute_forces:
        forces = atoms.get_forces()  # (N, 3) eV/Å

    symbols = list(atoms.get_chemical_symbols())

    return {
        "file": str(xyz_path),
        "n_atoms": len(atoms),
        "elements": symbols,
        "energy_eV": energy_eV,
        "energy_Ha": energy_eV * EV_TO_HARTREE,
        "energy_per_atom_eV": energy_eV / len(atoms),
        "forces_eV_Ang": forces,
    }


# ---------------------------------------------------------------------------
# Batch helper
# ---------------------------------------------------------------------------

def predict_batch(
    xyz_paths: list[str | Path],
    model_size: str = "small",
    compute_forces: bool = False,
    verbose: bool = True,
) -> list[dict]:
    """Predict energies for a list of .xyz files. Loads model once."""
    calc = load_model(model_size)
    results = []

    for path in xyz_paths:
        try:
            result = predict(path, calc, compute_forces=compute_forces)
            results.append(result)
            if verbose:
                _print_result(result)
        except Exception as exc:
            print(f"[WARN] {path}: {exc}", file=sys.stderr)

    return results


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _print_result(r: dict) -> None:
    unique_els = sorted(set(r["elements"]))
    print(
        f"\n{'-'*60}\n"
        f"File        : {Path(r['file']).name}\n"
        f"Composition : {' '.join(unique_els)}  ({r['n_atoms']} atoms)\n"
        f"Energy      : {r['energy_eV']:>14.6f} eV\n"
        f"            : {r['energy_Ha']:>14.6f} Hartree\n"
        f"Per atom    : {r['energy_per_atom_eV']:>14.6f} eV/atom"
    )
    if r["forces_eV_Ang"] is not None:
        f = r["forces_eV_Ang"]
        fmax = float(np.linalg.norm(f, axis=1).max())
        frms = float(np.sqrt((f**2).mean()))
        print(f"Max |force| : {fmax:.4f} eV/Ang  (RMS {frms:.4f} eV/Ang)")


def _print_summary(results: list[dict]) -> None:
    if len(results) < 2:
        return
    energies = np.array([r["energy_eV"] for r in results])
    print(
        f"\n{'='*60}\n"
        f"BATCH SUMMARY  ({len(results)} structures)\n"
        f"  Min energy : {energies.min():.6f} eV  ({Path(results[int(energies.argmin())]['file']).name})\n"
        f"  Max energy : {energies.max():.6f} eV\n"
        f"  Spread     : {energies.max() - energies.min():.6f} eV"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description="Predict DFT-quality energy from .xyz files using MACE-MP-0."
    )
    p.add_argument("xyz", nargs="+", help=".xyz file(s) to evaluate")
    p.add_argument(
        "--model",
        choices=["small", "medium", "large"],
        default="small",
        help="MACE-MP-0 model size (default: small ~30M params)",
    )
    p.add_argument(
        "--forces",
        action="store_true",
        help="Also compute atomic forces (eV/Å)",
    )
    p.add_argument(
        "--csv",
        metavar="FILE",
        help="Write results to a CSV file",
    )
    return p.parse_args()


def _write_csv(results: list[dict], path: str) -> None:
    import csv

    fieldnames = ["file", "n_atoms", "elements", "energy_eV", "energy_Ha", "energy_per_atom_eV"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            row = {k: r[k] for k in fieldnames}
            row["elements"] = " ".join(sorted(set(r["elements"])))
            w.writerow(row)
    print(f"\nResults saved to {path}")


def main():
    args = _parse_args()
    paths = [Path(p) for p in args.xyz]

    missing = [p for p in paths if not p.exists()]
    if missing:
        sys.exit(f"Files not found: {', '.join(str(m) for m in missing)}")

    print(f"MACE-MP-0 ({args.model})  |  {len(paths)} structure(s)")
    results = predict_batch(paths, model_size=args.model, compute_forces=args.forces)

    _print_summary(results)

    if args.csv:
        _write_csv(results, args.csv)


if __name__ == "__main__":
    main()
