"""Field specification for the experiment log: one source for templates, checks and schema.

A value is one of three things:

* blank (YAML `~` or nothing after the colon) — not recorded yet; `expt check` flags it;
* `n/a` — does not apply to this experiment;
* `unknown` — cannot be recovered;

or the value itself.  Blank is never read as zero, and neither sentinel is a gap.
"""
from __future__ import annotations

from dataclasses import dataclass, field

NA = "n/a"
UNKNOWN = "unknown"
SENTINELS = (NA, UNKNOWN)

STATUSES = ("planned", "ongoing", "done", "abandoned")
METHODS = ("solution", "ball_mill", "soxhlet", "recrystallization", "other")
ROLES = ("metal_source", "ligand", "additive", "solvent", "substrate")
KINDS = ("metal_salt", "ligand", "additive", "solvent", "other")


@dataclass(frozen=True)
class Leaf:
    """One value.  `kind`: num, int, temp, str, text, date, bool, path, enum, reagent, phase,
    entry."""

    kind: str
    help: str = ""
    unit: str = ""
    choices: tuple[str, ...] = ()
    required: bool = True


@dataclass(frozen=True)
class Block:
    fields: dict[str, Node]
    help: str = ""
    required: bool = True
    outcome: bool = False          # not asked for while the entry is `planned`


@dataclass(frozen=True)
class Items:
    item: Block
    help: str = ""
    required: bool = True
    outcome: bool = False
    key: str = ""                  # item field used to name an item in reports
    skeleton: tuple[dict, ...] = field(default=({},))   # blank items a new template carries


Node = Leaf | Block | Items


def _num(help: str = "", unit: str = "", required: bool = True) -> Leaf:
    return Leaf("num", help, unit, required=required)


def _temp(help: str = "") -> Leaf:
    return Leaf("temp", help, "°C, or RT")


def _str(help: str = "", required: bool = True) -> Leaf:
    return Leaf("str", help, required=required)


def _path(help: str = "data file, relative to the log folder (e.g. raw/2-01_pxrd.xy)",
          required: bool = True) -> Leaf:
    return Leaf("path", help, required=required)


REAGENT_LINE = Block({
    "code": Leaf("reagent", "code from reagents.yaml, or product:<id> for another entry's product"),
    "role": Leaf("enum", choices=ROLES),
    "mass_mg": _num("enter what was weighed; mmol is then computed", "mg", required=False),
    "volume_ml": _num("liquids and solvents", "mL", required=False),
    "mmol": _num("only when neither mass nor volume can give it", "mmol", required=False),
    "ratio": _num("planned molar equivalents (solvents: n/a)", "eq"),
    "note": _str(required=False),
})

CONDITIONS: dict[str, Block] = {
    "solution": Block({
        "temperature_c": _temp(),
        "time_h": _num("total reaction time", "h"),
        "addition": _str("what was added to what, how, and over how long"),
        "stirring": _str("e.g. magnetic, 500 rpm"),
        "atmosphere": _str("air, N2, Ar"),
    }),
    "ball_mill": Block({
        "mill": _str("model"),
        "frequency_hz": _num(unit="Hz"),
        "time_h": _num(unit="h"),
        "pauses": _str("e.g. 5 min every 30 min; none"),
        "jar_material": _str(),
        "jar_volume_ml": _num(unit="mL"),
        "ball_material": _str(),
        "ball_count": Leaf("int"),
        "ball_diameter_mm": _num(unit="mm"),
        "liquid": Leaf("reagent", "liquid-assisted grinding liquid (reagent code), or n/a"),
        "liquid_ul": _num("η is computed from this and the solid masses", "µL"),
        "atmosphere": _str("air, N2, Ar"),
    }),
    "soxhlet": Block({
        "thimble_contents": _str("reagent codes in the thimble"),
        "flask_contents": _str("reagent codes and solvent in the collection flask"),
        "solvent": _str("extraction solvent"),
        "temperature_c": _temp("bath or mantle temperature"),
        "time_h": _num(unit="h"),
        "cycles": Leaf("int", "siphon cycles, if counted", required=False),
    }),
    "recrystallization": Block({
        "solvent": _str("solvent system and ratio"),
        "solvent_ml": _num(unit="mL"),
        "temperature_c": _temp(),
        "time_h": _num(unit="h"),
        "cooling": _str("how it was cooled"),
        "vessel": _str("e.g. sealed vial, open beaker"),
    }),
    "other": Block({"description": Leaf("text", "the procedure")}),
}

WORKUP = Block({
    "washes": Items(Block({
        "solvent": _str(),
        "times": Leaf("int"),
        "mode": _str("e.g. vortex + centrifuge, filtration"),
    }), key="solvent"),
    "drying": Block({
        "method": _str("e.g. vacuum oven, hot plate, air"),
        "temperature_c": _temp(),
        "time_h": _num(unit="h"),
    }),
})

PRODUCT = Block({
    "appearance": _str("colour and form after drying"),
    "mass_mg": _num(unit="mg"),
    "yield_pct": _num("state the basis below", "%"),
    "yield_basis": _str("e.g. on Cu, assuming formula Cu2(ATF)"),
}, outcome=True)

CHARACTERIZATION = Block({
    "pxrd": Block({
        "phase": Leaf("phase", "label from phases.yaml"),
        "file": _path(),
        "instrument": _str("instrument, radiation, 2θ range"),
        "notes": _str("e.g. unreacted ligand peaks present", required=False),
    }),
    "ftir": Block({"file": _path(), "notes": _str(required=False)}),
    "chn": Block({
        "c_pct": _num(unit="wt%"),
        "h_pct": _num(unit="wt%"),
        "n_pct": _num(unit="wt%"),
        "file": _path(required=False),
    }),
    "tga": Block({
        "atmosphere": _str(),
        "file": _path(),
        "notes": _str("mass-loss steps", required=False),
    }),
    "bet": Block({"surface_m2_g": _num(unit="m²/g"), "file": _path(required=False)}),
    "microscopy": Block({
        "technique": _str("SEM, TEM, iDPC STEM"),
        "particle_size_nm": _num(unit="nm", required=False),
        "file": _path(required=False),
        "notes": _str(required=False),
    }),
    "oxidation_state": Block({
        "technique": _str("XPS, EPR"),
        "result": _str("e.g. Cu(II) only; Cu(I) present"),
        "file": _path(required=False),
    }),
    "other": _str("anything else measured", required=False),
}, outcome=True)

CELL = Block({
    "id": _str("e.g. c1; one item per cell built"),
    "date": Leaf("date", "YYYY-MM-DD"),
    "electrode": Block({
        "active_pct": _num(unit="wt%"),
        "carbon": _str("e.g. Super P"),
        "carbon_pct": _num(unit="wt%"),
        "binder": _str("e.g. PVDF"),
        "binder_pct": _num(unit="wt%"),
        "slurry_solvent": _str(),
        "current_collector": _str(),
        "loading_mg_cm2": _num("active material loading", "mg/cm²"),
        "area_cm2": _num(unit="cm²"),
    }),
    "cell": Block({
        "format": _str("e.g. CR2032"),
        "counter": _str("e.g. Li metal"),
        "electrolyte": _str("salt, concentration, solvent"),
        "electrolyte_ul": _num(unit="µL"),
        "separator": _str(),
        "temperature_c": _temp(),
    }),
    "protocol": Block({
        "v_min": _num(unit="V"),
        "v_max": _num(unit="V"),
        "rate_basis": _str("how current is set, e.g. mA/g of active material"),
        "formation": _str("formation cycles, or none"),
    }),
    "capacity_basis": Leaf("enum", "mass the capacity is divided by",
                           choices=("active", "electrode")),
    "carbon_blank": _str("cell or file giving the carbon+binder contribution, or n/a"),
    "results": Block({
        "capacity_mah_g": _num("reversible capacity", "mAh/g"),
        "capacity_rate_ma_g": _num("rate the capacity above was measured at", "mA/g"),
        "ice_pct": _num("first-cycle coulombic efficiency", "%"),
        "retention_pct": _num("capacity retained ...", "%"),
        "retention_cycles": Leaf("int", "... after this many cycles"),
        "retention_rate_ma_g": _num("... at this rate", "mA/g"),
        "ce_pct": _num("average coulombic efficiency over the cycling run", "%"),
        "rate_ladder": Items(Block({
            "rate_ma_g": _num(unit="mA/g"),
            "capacity_mah_g": _num(unit="mAh/g"),
        }), key="rate_ma_g"),
        "rate_recovery_pct": _num("capacity on returning to the first rate", "%"),
        "postmortem": _str("separator colour, dissolution, electrode state"),
    }),
    "files": Items(Block({"file": _path("cycler export, CV, EIS")}), key="file"),
})


def conditions_for(method: object) -> Block:
    return CONDITIONS.get(method if isinstance(method, str) else "", CONDITIONS["other"])


def entry_spec(method: object) -> Block:
    """The full entry for one synthesis method."""
    return Block({
        "id": Leaf("str", "matches the file name"),
        "title": _str("short description", required=False),
        "status": Leaf("enum", choices=STATUSES),
        "date": Leaf("date", "YYYY-MM-DD the synthesis started"),
        "operator": _str(),
        "notebook_ref": _str("lab notebook page"),
        "derived_from": Leaf("entry", "entry this one was copied from", required=False),
        "method": Leaf("enum", choices=METHODS),
        "reagents": Items(REAGENT_LINE, key="code", skeleton=(
            {"role": "metal_source"}, {"role": "ligand"}, {"role": "additive"},
            {"role": "solvent"})),
        "conditions": conditions_for(method),
        "workup": WORKUP,
        "product": PRODUCT,
        "characterization": CHARACTERIZATION,
        "cells": Items(CELL, outcome=True, key="id"),
        "notes": Leaf("text", required=False),
    })


# Reagents and phases are keyed by code / label; these are the per-item fields.
REAGENT = Block({
    "name": _str("full name"),
    "kind": Leaf("enum", choices=KINDS),
    "formula": _str("with hydrate, e.g. Cu(C2H3O2)2·H2O; molar mass is computed from it"),
    "mw": _num("only if the formula cannot give it", "g/mol", required=False),
    "cas": _str(),
    "smiles": _str("ligands and additives at least", required=False),
    "supplier": _str(),
    "purity": _str(),
    "density_g_ml": _num("liquids", "g/mL", required=False),
    "concentration": _str("solutions, e.g. 28 wt% NH3", required=False),
})

PHASE = Block({
    "description": _str(),
    "peaks_2theta": _str("defining reflections, e.g. 6.2, 12.4, 26.8", required=True),
    "radiation": _str("e.g. Cu Kα"),
    "reference_file": _path(required=False),
})

OUTCOME_SECTIONS = tuple(k for k, v in entry_spec("other").fields.items()
                         if getattr(v, "outcome", False))
