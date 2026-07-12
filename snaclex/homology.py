"""Homology-based structure selection (Stage 4, tier 2).

When a protein has no *direct* experimental structure, structurally-characterized
homologs are the next-best evidence — above a predicted model in the hierarchy
(experimental → experimental-homologous → predicted → sequence-only). This module
turns RCSB sequence-search hits into homologous-structure records, each carrying
its sequence identity and an explicit "not the requested protein" caveat, and
gates docking so a distant homolog is not treated as a reliable receptor.

The parsing/classification is pure (``classify_hits``); only ``rcsb.sequence_search``
touches the network, so this is offline-testable by feeding recorded hits.
"""

from __future__ import annotations

import datetime

from . import rcsb
from .http_util import FetchError

# A homolog needs to be reasonably close before we'd even *offer* it for docking.
_DOCKING_IDENTITY_MIN = 0.40
# At/above this we flag the hit as possibly the same protein (near-identical),
# not a distinct homolog, so identity context is never silently lost.
_NEAR_IDENTICAL = 0.95


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def classify_hits(hits: list[dict], *, direct_pdb_ids=(), max_report: int = 10) -> list[dict]:
    """Turn raw sequence-search hits into homologous-structure records (pure).

    Excludes the protein's own (direct) PDB entries; keeps 1:N (a homolog may
    appear as several entities) by de-duplicating on PDB id. Identity context is
    always preserved — never collapsed.
    """
    direct = {p.upper() for p in direct_pdb_ids}
    seen: set[str] = set()
    out: list[dict] = []
    for h in hits:
        pid = (h.get("pdb_id") or "").upper()
        if not pid or pid in direct or pid in seen:
            continue
        seen.add(pid)
        ident = h.get("identity")
        near_identical = ident is not None and ident >= _NEAR_IDENTICAL
        dockable = ident is not None and ident >= _DOCKING_IDENTITY_MIN
        out.append({
            "kind": "homologous", "source": "RCSB", "pdb_id": pid,
            "entity": h.get("entity"), "sequence_identity": ident,
            "evalue": h.get("evalue"),
            "near_identical": near_identical,
            "docking_suitable": bool(dockable),
            "docking_caveat": (
                "homolog of the requested protein, not the protein itself — "
                "binding-site residues may differ; transfer with caution"),
            "provenance": {
                "source": "RCSB PDB (sequence search)", "source_id": pid,
                "retrieved_utc": _now(), "evidence_category": "experimental",
                "limitations": "structural homolog; sequence identity < 100%"
                               if not near_identical else "near-identical sequence match",
            },
        })
        if len(out) >= max_report:
            break
    return out


def find_homologous_structures(sequence: str, *, direct_pdb_ids=(),
                               identity_cutoff: float = 0.3,
                               limit: int = 25, max_report: int = 10) -> dict:
    """Search RCSB for structural homologs of a sequence and classify them.

    Never raises for a routine miss/short sequence; returns ``available`` so the
    caller can render an explicit state.
    """
    seq = "".join((sequence or "").split())
    if len(seq) < 20:
        return {"available": False, "reason": "sequence too short for a structural search",
                "structures": []}
    try:
        hits = rcsb.sequence_search(seq, identity_cutoff=identity_cutoff, limit=limit)
    except FetchError as exc:
        return {"available": False, "reason": str(exc), "structures": []}
    structures = classify_hits(hits, direct_pdb_ids=direct_pdb_ids, max_report=max_report)
    return {
        "available": True,
        "structures": structures,
        "method": {
            "service": "RCSB sequence search (mmseqs2)",
            "identity_cutoff": identity_cutoff,
            "docking_identity_min": _DOCKING_IDENTITY_MIN,
            "note": "Homologous experimental structures rank above predicted "
                    "models but below a direct experimental structure.",
        },
    }
