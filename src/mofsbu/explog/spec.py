"""Field specification for the experiment log: one source for the form, the check and the views.

A value is one of three things:

* blank (absent or null) — not recorded yet; the check flags it;
* `n/a` — does not apply to this experiment;
* `unknown` — cannot be recovered;

or the value itself.  Blank is never read as zero, and neither sentinel is a gap.

A recipe is a list of *preparations* (each solution, suspension or solid mixture made
before the reaction: what is in it, in which solvents, and how it was mixed), then an
ordered list of *steps* (what was combined or done, how, for how long, at what
temperature, with what mixing), then the *apparatus* the method needs.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

NA = "n/a"
UNKNOWN = "unknown"
SENTINELS = (NA, UNKNOWN)

STATUSES = ("planned", "ongoing", "done", "abandoned")
METHODS = ("solution", "solvothermal", "ball_mill", "soxhlet", "recrystallization", "other")
ROLES = ("metal_source", "ligand", "additive", "solvent", "substrate")
KINDS = ("metal_salt", "ligand", "additive", "solvent", "other")
FORMS = ("solution", "suspension", "solid_mixture", "neat")
ACTIONS = ("add", "combine", "stir", "sonicate", "heat", "cool", "reflux", "mill", "extract",
           "age", "other")


@dataclass(frozen=True)
class Leaf:
    """One value.  `kind`: num, int, temp, str, text, date, enum, path, reagent, phase, entry."""

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
    outcome: bool = False          # not asked for while the entry is `planned` or `abandoned`


@dataclass(frozen=True)
class Items:
    item: Block
    help: str = ""
    required: bool = True
    outcome: bool = False
    key: str = ""                  # item field used to name an item in reports


Node = Leaf | Block | Items


def _num(help: str = "", unit: str = "", required: bool = True) -> Leaf:
    return Leaf("num", help, unit, required=required)


def _temp(help: str = "", required: bool = True) -> Leaf:
    return Leaf("temp", help, "°C, or RT", required=required)


def _str(help: str = "", required: bool = True) -> Leaf:
    return Leaf("str", help, required=required)


def _path(help: str = "data file, relative to the log folder (e.g. raw/2-01_pxrd.xy)",
          required: bool = True) -> Leaf:
    return Leaf("path", help, required=required)


COMPONENT = Block({
    "code": Leaf("reagent", "reagent code, or product:<id> for another entry's product"),
    "role": Leaf("enum", choices=ROLES),
    "mass_mg": _num("what was weighed; mmol is computed from it", "mg", required=False),
    "volume_ml": _num("liquids, solvents and stock solutions", "mL", required=False),
    "mmol": _num("only when neither mass nor volume can give it", "mmol", required=False),
    "ratio": _num("planned molar equivalents (not asked for solvents)", "eq"),
    "note": _str("e.g. few drops; added last", required=False),
})

PREPARATION = Block({
    "label": _str("e.g. A: salt solution"),
    "form": Leaf("enum", choices=FORMS),
    "components": Items(COMPONENT, "solutes and solvents in this preparation", key="code"),
    "mix_method": _str("how it was made up: magnetic stirring, sonication, vortex, grinding, "
                       "none"),
    "mix_speed_rpm": _num(unit="rpm", required=False),
    "mix_time_min": _num(unit="min"),
    "temperature_c": _temp(),
    "observation": _str("e.g. fully dissolved; cloudy", required=False),
})

STEP = Block({
    "action": Leaf("enum", choices=ACTIONS),
    "what": _str("e.g. B into A; jar sealed and milled"),
    "addition": _str("dropwise, portionwise, at once — and over how long", required=False),
    "duration_h": _num("how long this step lasted", "h"),
    "temperature_c": _temp(),
    "mixing": _str("magnetic stirring, overhead, shaking, none"),
    "speed_rpm": _num(unit="rpm", required=False),
    "atmosphere": _str("air, N2, Ar"),
    "note": _str(required=False),
})

APPARATUS: dict[str, Block] = {
    "solution": Block({
        "vessel": _str("e.g. 250 mL round-bottom flask"),
    }),
    "solvothermal": Block({
        "vessel": _str("e.g. 23 mL Teflon-lined autoclave"),
        "fill_pct": _num(unit="%"),
        "oven": _str(),
        "heating_rate": _str("e.g. 2 °C/min, or placed in a hot oven"),
        "cooling": _str("e.g. 0.1 °C/min, or natural"),
    }),
    "ball_mill": Block({
        "mill": _str("model"),
        "frequency_hz": _num(unit="Hz"),
        "jar_material": _str(),
        "jar_volume_ml": _num(unit="mL"),
        "ball_material": _str(),
        "ball_count": Leaf("int"),
        "ball_diameter_mm": _num(unit="mm"),
        "liquid": Leaf("reagent", "liquid-assisted grinding liquid, or n/a"),
        "liquid_ul": _num("η is computed from this and the solid masses", "µL"),
        "pauses": _str("e.g. 5 min every 30 min; none"),
    }),
    "soxhlet": Block({
        "thimble_contents": _str("what sat in the thimble"),
        "flask_contents": _str("what was in the collection flask"),
        "thimble": _str("thimble material and size", required=False),
        "cycles": Leaf("int", "siphon cycles, if counted", required=False),
    }),
    "recrystallization": Block({
        "vessel": _str("e.g. sealed 20 mL vial"),
        "cooling": _str("how it was cooled"),
    }),
    "other": Block({"description": Leaf("text", "the apparatus")}),
}

WORKUP = Block({
    "separation": _str("e.g. centrifuge 8000 rpm 10 min; filtration"),
    "washes": Items(Block({
        "solvent": Leaf("reagent"),
        "volume_ml": _num("per wash", "mL"),
        "times": Leaf("int"),
        "mode": _str("e.g. vortex, soak 1 h", required=False),
    }), key="solvent"),
    "drying": Block({
        "method": _str("e.g. vacuum oven, hot plate, air"),
        "temperature_c": _temp(),
        "time_h": _num(unit="h"),
    }),
    "storage": _str("e.g. desiccator, glovebox", required=False),
})

PRODUCT = Block({
    "appearance": _str("colour and form after drying"),
    "mass_mg": _num(unit="mg"),
    "yield_pct": _num("state the basis below", "%"),
    "yield_basis": _str("e.g. on Cu, assuming formula Cu2(ATF)"),
}, outcome=True)

CHARACTERIZATION = Block({
    "pxrd": Block({
        "phase": Leaf("phase", "label from the phase table"),
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


def apparatus_for(method: object) -> Block:
    return APPARATUS.get(method if isinstance(method, str) else "", APPARATUS["other"])


def entry_spec(method: object) -> Block:
    """The full entry for one synthesis method."""
    return Block({
        "id": Leaf("str", "letters, digits, . _ - ; it names the file"),
        "title": _str("short description", required=False),
        "status": Leaf("enum", choices=STATUSES),
        "date": Leaf("date", "the day the synthesis started"),
        "operator": _str(),
        "notebook_ref": _str("lab notebook page"),
        "derived_from": Leaf("entry", "entry this one was copied from", required=False),
        "method": Leaf("enum", choices=METHODS),
        "preparations": Items(PREPARATION, "each solution, suspension or solid mixture",
                              key="label"),
        "steps": Items(STEP, "in the order they were done", key="action"),
        "apparatus": apparatus_for(method),
        "workup": WORKUP,
        "product": PRODUCT,
        "characterization": CHARACTERIZATION,
        "cells": Items(CELL, "one per coin cell built", outcome=True, key="id"),
        "notes": Leaf("text", "anything the fields do not hold", required=False),
    })


# Reagents and phases are tables keyed by code / label; these are the per-row fields.
REAGENT = Block({
    "name": _str("full name"),
    "kind": Leaf("enum", choices=KINDS),
    "formula": _str("with hydrate, e.g. Cu(C2H3O2)2·H2O; molar mass is computed from it"),
    "mw": _num("only if the formula cannot give it", "g/mol", required=False),
    "cas": _str(),
    "smiles": _str("ligands and additives at least", required=False),
    "supplier": _str(),
    "purity": _str(),
    "density_g_ml": _num("liquids and solutions; turns a volume into a mass", "g/mL",
                         required=False),
    "wt_pct": _num("solutions sold by assay, e.g. 28 for 28 wt% NH3", "wt%", required=False),
    "molarity_m": _num("stock solutions made up to a concentration", "mol/L", required=False),
})

PHASE = Block({
    "description": _str(),
    "peaks_2theta": _str("defining reflections, e.g. 6.2, 12.4, 26.8"),
    "radiation": _str("e.g. Cu Kα"),
    "reference_file": _path(required=False),
})

OUTCOME_SECTIONS = tuple(k for k, v in entry_spec("other").fields.items()
                         if getattr(v, "outcome", False))


def describe(node: Node) -> dict:
    """The spec as plain data, for the page that builds its form from it."""
    if isinstance(node, Leaf):
        return {"type": "leaf", **asdict(node), "choices": list(node.choices)}
    if isinstance(node, Block):
        return {"type": "block", "help": node.help, "required": node.required,
                "outcome": node.outcome,
                "fields": [[k, describe(v)] for k, v in node.fields.items()]}
    return {"type": "items", "help": node.help, "required": node.required,
            "outcome": node.outcome, "key": node.key, "item": describe(node.item)}


def describe_all() -> dict:
    return {
        "sentinels": list(SENTINELS),
        "entry": describe(entry_spec("other")),
        "apparatus": {m: describe(b) for m, b in APPARATUS.items()},
        "reagent": describe(REAGENT),
        "phase": describe(PHASE),
    }
