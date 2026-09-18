"""The choice vector: the discrete record of how a structure was built (D13, M5/S1).

D13 says `construct` is a deterministic function of `(choice_vector, seed)` and **emits its
choice vector**; D11 makes that vector the L3 label — conformer identity is
provenance-primary and geometry only verifies.  Both claims rest on one operation: turning
the vector into a short, stable key that two machines agree on.  That is this module.

It is deliberately not a new vocabulary.  `geometry.placer.place_mononuclear` has been
emitting the shape since rev 18 — ``{metal, geometry, cn, d_ml, seed, ligands: [{ligand,
mode, torsion_well, azimuth_step, oop_step, ...}], donor_distances}`` — and
`azimuth_step` / `oop_step` are already "a stochastic search reduced to a recorded index",
which is exactly what D13 asks a conformer coordinate to be.  So a `ChoiceVector` wraps a
plain mapping rather than replacing it, and every producer stays a producer of dicts.

**What the digest covers, and what it deliberately does not:**

* **The seed is out.**  The schema gives `seed` its own column next to
  `choice_vector_digest` for a reason: a seed is a Kind-A coordinate (§6.7) — stochastic
  noise to be deduplicated — not a Kind-B branch.  Two geometries built from one choice
  vector at different seeds are *the same L3 branch, sampled twice*, and the query
  `ix_geometries_choice` exists to serve ("give me every geometry from this choice") is
  precisely the one that must return both so clustering can collapse them.  Putting the
  seed in the key would make that family unfindable and every duplicate look like a
  distinct conformer.
* **A `digest` key is out**, because a vector that carries its own digest cannot be
  digested twice to the same value, and stored JSON round-trips through exactly that path.
* **Floats are quantised** to `FLOAT_DP` decimals.  The values here are angles and
  distances in degrees and angstrom — the placer already rounds most of them to 3 dp — so
  1e-6 is far below anything chemically or numerically meaningful, while float arithmetic
  that differs in the last bits between two machines is a real and boring way to lose
  reproducibility (ground rule 9).

**Why the version is a prefix and not a column.**  The L1 key keeps its recipe version in a
separate column because `block_id` is an address and an address may not carry incidental
text.  A choice digest is not part of any address, and there is no column for it — so it
carries `ALGO_VERSIONS["choice_vector"]` in front, which makes a stale digest
self-identifying and makes it impossible for a future recipe to produce a value that
silently compares equal to this one's.

Bump that version when anything below this line changes — and also when a PRODUCER starts
recording a coordinate it did not record before, which is the case that looks like it is
somebody else's business and is not: a key over a larger set of choices is a different
key, whatever the canonicalisation did.  Old rows keep their prefix and are never
re-labelled (ground rule 6).
"""
from __future__ import annotations

import hashlib
import json
import numbers
from dataclasses import dataclass, field
from typing import Any, Mapping

from mofsbu._types import MofsbuError
from mofsbu.versions import ALGO_VERSIONS

#: Decimals kept when a float enters the key.  Part of the recipe: changing it changes
#: every digest, so it moves with `ALGO_VERSIONS["choice_vector"]`.
FLOAT_DP = 6

#: Top-level keys dropped before hashing.  Not nested ones — a `seed` inside a per-ligand
#: record would mean something else (nothing emits one today), and silently eating keys by
#: name at every depth is how a real choice disappears from a key.
EXCLUDED = ("seed", "digest")


class ChoiceVectorError(MofsbuError):
    """A choice vector that cannot be turned into a stable key."""


def _canonical_value(value: Any, path: str) -> Any:
    """One value, in the form the digest sees.  Raises rather than guessing."""
    if value is None or isinstance(value, str):
        return value
    if isinstance(value, bool):
        # Before the Integral check: `bool` IS an `Integral`, and coercing True to 1
        # would make `{"estimated": True}` and `{"estimated": 1}` the same key while
        # they round-trip through JSON as different types.
        return value
    if isinstance(value, numbers.Integral):
        return int(value)                     # numpy int64 is Integral but not int
    if isinstance(value, numbers.Real):
        f = float(value)
        if f != f or f in (float("inf"), float("-inf")):
            raise ChoiceVectorError(
                f"{path} is {value!r}. A non-finite number in a choice vector means the "
                "search that produced it failed; hashing it would store that failure as "
                "a reproducible choice.")
        r = round(f, FLOAT_DP)
        return r + 0.0                        # -0.0 and 0.0 are one value here
    if isinstance(value, Mapping):
        return {str(k): _canonical_value(v, f"{path}.{k}") for k, v in sorted(
            value.items(), key=lambda kv: str(kv[0]))}
    if isinstance(value, (list, tuple)):
        # Order is meaningful — ligand 0 and ligand 1 are different choices — so lists
        # are NOT sorted.  Only the tuple/list distinction is flattened, because JSON
        # has one sequence type and a stored vector must digest the same after a
        # round-trip through it.
        return [_canonical_value(v, f"{path}[{i}]") for i, v in enumerate(value)]
    raise ChoiceVectorError(
        f"{path} is a {type(value).__name__}, which has no JSON form. A choice vector is "
        "stored and replayed as JSON, so every value in it must survive that trip.")


def canonical(data: Mapping[str, Any]) -> dict[str, Any]:
    """The comparable form of a choice vector: sorted, quantised, JSON-safe, seed-free."""
    if not isinstance(data, Mapping):
        raise ChoiceVectorError(
            f"a choice vector is a mapping of named choices, not {type(data).__name__}")
    return {str(k): _canonical_value(v, str(k))
            for k, v in sorted(data.items(), key=lambda kv: str(kv[0]))
            if str(k) not in EXCLUDED}


def canonical_json(data: Mapping[str, Any]) -> str:
    """The exact bytes the digest is taken over.  This is what gets stored."""
    return json.dumps(canonical(data), sort_keys=True, separators=(",", ":"),
                      allow_nan=False)


def digest_of(data: Mapping[str, Any]) -> str:
    """`<version>:<sha256>` over the canonical form — what `choice_vector_digest` holds."""
    payload = canonical_json(data).encode()
    return f"{ALGO_VERSIONS['choice_vector']}:{hashlib.sha256(payload).hexdigest()}"


@dataclass(frozen=True, eq=False)
class ChoiceVector:
    """One path through the construction tree (§6.7), as a record with a stable key.

    Wraps a mapping instead of naming fields, because the key set is open by design: a
    mononuclear placement's choices are not a join's are not a protomer branch's, and
    §6.7 expects all three to be leaves of the same tree.  What is fixed is the *recipe* —
    what gets dropped, how floats are quantised, how it serialises — and that lives here
    so every producer inherits it rather than re-deciding.

    `digest` is recomputed on every read rather than cached at construction: `data` is a
    plain dict and Python cannot stop a caller mutating it, so a cached key would go
    quietly stale exactly where staleness is unrecoverable.  The vectors are small.
    """

    data: dict[str, Any] = field(default_factory=dict)
    #: The producer's seed.  Kept OUT of the key (see the module docstring) but not
    #: thrown away: the placer has always written one into its dict, and a caller that
    #: stores a geometry without passing `seed=` separately would otherwise lose the one
    #: number `construct` needs to reproduce the coordinates.
    seed: int | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        raw = self.data
        object.__setattr__(self, "data", canonical(raw))
        if self.seed is None:
            object.__setattr__(self, "seed", _seed_in(raw))

    # -- identity ----------------------------------------------------------

    @property
    def digest(self) -> str:
        return digest_of(self.data)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, ChoiceVector):
            return self.data == other.data
        if isinstance(other, Mapping):
            return self.data == canonical(other)
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.digest)

    def __repr__(self) -> str:
        return f"ChoiceVector({self.digest[:16]}…, {len(self.data)} choices)"

    # -- storage -----------------------------------------------------------

    def to_json(self) -> str:
        return canonical_json(self.data)

    @classmethod
    def from_json(cls, text: str) -> "ChoiceVector":
        """Rebuild a vector from a stored `choice_vector_json`.

        This is the half of "replay" that belongs here: handing back the vector that
        `construct` is a deterministic function of.  Turning it back into coordinates is
        `assembly.construct`'s job (M5/S4), and needs the row's `seed` as well — which is
        why the seed is a column rather than a key inside the vector.
        """
        loaded = json.loads(text)
        if not isinstance(loaded, dict):
            raise ChoiceVectorError(
                f"stored choice vector is a {type(loaded).__name__}, not an object")
        return cls(loaded)

    @classmethod
    def coerce(cls, value: "ChoiceVector | Mapping[str, Any] | None") -> "ChoiceVector | None":
        """Accept whichever form a caller has; `None` stays `None`.

        Producers emit plain dicts (the placer does) and the registry should not care.
        """
        if value is None:
            return None
        return value if isinstance(value, cls) else cls(value)


def _seed_in(data: Any) -> int | None:
    """The seed a producer wrote into its vector, if it wrote a usable one."""
    if not isinstance(data, Mapping):
        return None
    value = data.get("seed")
    if value is None or isinstance(value, bool) or not isinstance(value, numbers.Integral):
        return None
    return int(value)
