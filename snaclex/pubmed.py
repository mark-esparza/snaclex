"""PubMed literature adapter (optional enrichment, Phase 3).

The record's Literature tab is populated for free from UniProt's *curated*
references (see ``uniprot._references``). This adapter adds optional deeper
enrichment via NCBI E-utilities — ``elink`` to find PubMed articles linked to a
protein accession, and ``esummary`` for article metadata — reusing the same
E-utilities plumbing and ``NCBI_API_KEY`` rate policy as ``ncbi.py``.

It is opt-in per request and degrades gracefully; a failure never breaks the core
record. Parsing is pure and offline-testable. Only public citation *metadata* is
retrieved (no full text), so there is no redistribution concern.
"""

from __future__ import annotations

import os
import urllib.parse

from .http_util import FetchError, fetch_json

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _api_key_suffix() -> str:
    key = os.environ.get("NCBI_API_KEY")
    return f"&api_key={urllib.parse.quote(key)}" if key else ""


def elink_pmids(accession: str, *, limit: int = 20) -> list[str]:
    """PubMed ids linked to a protein accession (elink protein→pubmed)."""
    acc = urllib.parse.quote((accession or "").strip())
    if not acc:
        return []
    url = (f"{_EUTILS}/elink.fcgi?dbfrom=protein&db=pubmed&id={acc}"
           f"&retmode=json{_api_key_suffix()}")
    try:
        data = fetch_json(url)
    except FetchError:
        return []
    return parse_elink(data)[:limit]


def fetch_summaries(pmids: list[str]) -> list[dict]:
    """Article metadata for a list of PMIDs (esummary)."""
    ids = [str(p) for p in pmids if p]
    if not ids:
        return []
    url = (f"{_EUTILS}/esummary.fcgi?db=pubmed&id={urllib.parse.quote(','.join(ids))}"
           f"&retmode=json{_api_key_suffix()}")
    try:
        data = fetch_json(url)
    except FetchError:
        return []
    return parse_esummary(data)


def search(term: str, *, limit: int = 20) -> list[str]:
    """Topical PubMed search → PMIDs (esearch)."""
    q = urllib.parse.quote((term or "").strip())
    if not q:
        return []
    url = (f"{_EUTILS}/esearch.fcgi?db=pubmed&term={q}&retmax={int(limit)}"
           f"&retmode=json{_api_key_suffix()}")
    try:
        data = fetch_json(url)
    except FetchError:
        return []
    return ((data.get("esearchresult") or {}).get("idlist")) or []


# ---------------------------------------------------------------------------
# Pure parsing
# ---------------------------------------------------------------------------
def parse_elink(data: dict) -> list[str]:
    """Extract linked PMIDs from an elink JSON response."""
    pmids: list[str] = []
    for linkset in (data or {}).get("linksets") or []:
        for db in linkset.get("linksetdbs") or []:
            for link in db.get("links") or []:
                sid = str(link)
                if sid not in pmids:
                    pmids.append(sid)
    return pmids


def parse_esummary(data: dict) -> list[dict]:
    """Normalize a pubmed esummary payload into citation dicts.

    VERIFY: field names (``source``, ``pubdate``, ``authors[].name``,
    ``articleids[].idtype == 'doi'``) follow the documented esummary shape.
    """
    result = (data or {}).get("result") or {}
    out: list[dict] = []
    for uid in result.get("uids") or []:
        rec = result.get(uid) or {}
        authors = [a.get("name") for a in rec.get("authors") or [] if a.get("name")]
        doi = None
        for aid in rec.get("articleids") or []:
            if aid.get("idtype") == "doi":
                doi = aid.get("value")
        out.append({
            "pmid": uid,
            "title": rec.get("title"),
            "authors": ", ".join(authors[:3]) + (" et al." if len(authors) > 3 else ""),
            "journal": rec.get("source"),
            "year": (rec.get("pubdate") or "")[:4] or None,
            "doi": doi,
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/",
            "source": "PubMed (elink)",
        })
    return out
