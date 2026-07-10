"""InterPro adapter — families, domains, repeats and important sites.

An *optional, independent* enrichment adapter (Phase 2): it is env-gated
(``SNACLEX_ENABLE_INTERPRO``) and degrades gracefully, so a failure, a rate
limit, or it simply being disabled never breaks the core sequence-first record —
matching the repo's "opt-in enrichment track" philosophy (ROADMAP.md).

Parsing (``parse_entries``) is pure and offline-testable. Field paths follow the
documented InterPro API shape and are marked VERIFY until round-tripped against a
live response (blocked by egress policy in CI). License: InterPro data are
CC BY 4.0 — attribute InterPro / EMBL-EBI.
"""

from __future__ import annotations

import datetime
import os
import urllib.parse

from .http_util import FetchError, fetch_json

_BASE = "https://www.ebi.ac.uk/interpro/api"


def enabled() -> bool:
    return os.environ.get("SNACLEX_ENABLE_INTERPRO", "").lower() in (
        "1", "true", "on", "yes")


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fetch_domains(uniprot_accession: str) -> dict:
    """Return InterPro entries mapped onto a UniProt accession (or 'unavailable').

    Never raises for a routine miss/disabled state — returns a dict with
    ``available`` so the caller can render an explicit state, not a blank.
    """
    if not enabled():
        return {"available": False,
                "reason": "InterPro enrichment disabled (set SNACLEX_ENABLE_INTERPRO=1)",
                "entries": []}
    acc = urllib.parse.quote((uniprot_accession or "").strip())
    if not acc:
        return {"available": False, "reason": "empty accession", "entries": []}
    url = f"{_BASE}/entry/all/protein/uniprot/{acc}/?page_size=100"
    try:
        data = fetch_json(url)
    except FetchError as exc:
        return {"available": False, "reason": str(exc), "entries": []}
    return {
        "available": True,
        "entries": parse_entries(data),
        "provenance": {
            "source": "InterPro", "source_id": uniprot_accession,
            "source_version": None, "retrieved_utc": _now(),
            "evidence_category": "curated",
            "limitations": "InterPro member-database integration; CC BY 4.0",
        },
    }


def parse_entries(data: dict) -> list[dict]:
    """Normalize the InterPro protein-entry response into domain records.

    VERIFY: paths ``results[].metadata.{accession,name,type,source_database}`` and
    ``results[].proteins[].entry_protein_locations[].fragments[].{start,end}``.
    """
    out: list[dict] = []
    for r in (data or {}).get("results", []):
        md = r.get("metadata") or {}
        locations: list[dict] = []
        for prot in r.get("proteins") or []:
            for epl in prot.get("entry_protein_locations") or []:
                for frag in epl.get("fragments") or []:
                    if frag.get("start") is not None:
                        locations.append({"start": frag.get("start"),
                                          "end": frag.get("end")})
        out.append({
            "accession": md.get("accession"),
            "name": md.get("name"),
            "type": (md.get("type") or "").lower() or None,
            "source_database": md.get("source_database"),
            "locations": locations,
            "source": "InterPro",
        })
    return out
