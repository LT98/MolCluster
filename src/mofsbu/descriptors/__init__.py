"""descriptors: tabulated per-donor and per-metal descriptors.

pKa, HSAB class, formal charge, ionic radius, exchange lability — the shared substrate
for C5/C6/C7 (DESIGN §6.6).  `tables` is built and loaded; `ease` is the policy layer on
top of it and is still a stub, because C5 and C7 are open checkpoints.
"""
from __future__ import annotations

from mofsbu.descriptors.ease import EaseRecord, activation_ease, hsab_match
from mofsbu.descriptors.tables import (
    DescriptorTableError, DonorDescriptor, MetalDescriptor, UnknownDescriptor, donor,
    donor_table, known_donor_types, known_ions, metal, metal_table, sync_to_registry,
)

__all__ = [
    "DescriptorTableError", "DonorDescriptor", "EaseRecord", "MetalDescriptor",
    "UnknownDescriptor", "activation_ease", "donor", "donor_table", "hsab_match",
    "known_donor_types", "known_ions", "metal", "metal_table", "sync_to_registry",
]
