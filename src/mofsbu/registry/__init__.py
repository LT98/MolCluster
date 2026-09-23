"""SQLite registry + content-addressed blob store."""
from __future__ import annotations

from mofsbu.registry.api import (
    MethodSpec, Provenance, Put, RegistryError, alias_fragment, backfill_choice_digests,
    best_geometry, canonical_map,
    catalog_drift,
    display_label,
    export_xyz, find, formula, geometries_from_choice, geometry_xyz, get_graph,
    get_site_state, get_sites,
    get_structure, find_method_id, incoming_routes, method_id,
    fragment_aliases, put_geometry, put_reaction, put_site_state, put_sites, put_structure,
    put_solvation_correction, relabel_all, set_hidden, solvation_corrections,
)
from mofsbu.registry.db import ReadOnlyRegistry, Registry, ensure_registry
from mofsbu.registry.store import BlobStore

__all__ = [
    "BlobStore", "MethodSpec", "Provenance", "Put", "ReadOnlyRegistry", "Registry",
    "RegistryError",
    "backfill_choice_digests", "best_geometry", "canonical_map", "catalog_drift",
    "display_label", "ensure_registry", "export_xyz", "find", "formula",
    "find_method_id", "geometries_from_choice", "geometry_xyz", "get_graph",
    "get_site_state", "get_sites",
    "get_structure", "incoming_routes", "method_id", "put_geometry",
    "put_reaction", "put_site_state", "put_sites", "put_structure", "relabel_all",
    "alias_fragment", "fragment_aliases", "set_hidden",
    "put_solvation_correction", "solvation_corrections",
]
