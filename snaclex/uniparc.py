"""UniParc adapter — sequence archive & provenance layer.

UniParc is used *only* for sequence-level identity, deduplication, cross-database
tracking, historical versions, and detecting obsolete/replaced accessions — never
as a functional-annotation source (see task + docs/platform/03).

The key bridge: a raw sequence's CRC-64 (``snaclex.checksum.crc64``) is queried
against UniParc ``checksum:`` to recover the stable UPI and every database record
that shares that exact sequence.
"""

from __future__ import annotations

import datetime
import urllib.parse

from .http_util import FetchError, fetch_json

_BASE = "https://rest.uniprot.org/uniparc"


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _provenance(upi) -> dict:
    return {
        "source": "UniParc",
        "source_id": upi,
        "source_version": None,
        "retrieved_utc": _now(),
        "evidence_category": "curated",
        "limitations": "sequence archive only — no functional annotation",
    }


# ---------------------------------------------------------------------------
# Network entry points
# ---------------------------------------------------------------------------
def fetch_entry(upi: str) -> dict:
    upi = urllib.parse.quote((upi or "").strip())
    if not upi:
        raise FetchError("empty UPI")
    return parse_entry(fetch_json(f"{_BASE}/{upi}.json"))


def search_by_checksum(crc64: str, *, limit: int = 5) -> list[dict]:
    """Find UniParc entries whose sequence CRC-64 matches (exact-sequence bridge)."""
    url = (f"{_BASE}/search?query=checksum:{urllib.parse.quote(crc64)}"
           f"&format=json&size={int(limit)}")
    data = fetch_json(url)
    return [parse_entry(r)["data"] for r in (data or {}).get("results", [])]


def search_by_uniprot(accession: str, *, limit: int = 5) -> list[dict]:
    url = (f"{_BASE}/search?query=uniprot:{urllib.parse.quote(accession)}"
           f"&format=json&size={int(limit)}")
    data = fetch_json(url)
    return [parse_entry(r)["data"] for r in (data or {}).get("results", [])]


# ---------------------------------------------------------------------------
# Pure parsing
# ---------------------------------------------------------------------------
def parse_entry(data: dict) -> dict:
    upi = data.get("uniParcId")
    seq = data.get("sequence") or {}
    xrefs = []
    for x in data.get("uniParcCrossReferences") or []:
        xrefs.append({
            "database": x.get("database"),
            "id": x.get("id"),
            "active": x.get("active"),
            "version": x.get("version"),
            "version_seq": x.get("versionI"),
            "created": x.get("created"),
            "last_updated": x.get("lastUpdated"),
        })
    active_dbs = sorted({x["database"] for x in xrefs if x.get("active") and x.get("database")})
    obsolete = [x for x in xrefs if x.get("active") is False]
    fragment = {
        "upi": upi,
        "sequence": {
            "value": seq.get("value"),
            "length": seq.get("length"),
            "crc64": seq.get("crc64"),
            "md5": seq.get("md5"),
        },
        "cross_references": xrefs,
        "active_databases": active_dbs,
        "obsolete_references": obsolete,
    }
    return {"data": fragment, "provenance": _provenance(upi), "status": "ok"}


def uniprot_accessions(fragment: dict) -> list[str]:
    """Pull UniProtKB/TrEMBL accessions out of a parsed UniParc fragment."""
    out: list[str] = []
    for x in fragment.get("cross_references") or []:
        db = (x.get("database") or "").upper()
        if db.startswith("UNIPROT") and x.get("id") and x["id"] not in out:
            out.append(x["id"])
    return out
