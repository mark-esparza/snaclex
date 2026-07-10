"""NCBI Protein / RefSeq adapter via E-utilities.

Used for protein sequence records, reference proteins, organism/strain, and
gene/nucleotide/PubMed links. Honors NCBI's rate policy: anonymous 3 req/s,
10 req/s when ``NCBI_API_KEY`` is set (the shared HTTP layer adds retry/backoff;
the key is appended here). RefSeq accession classes (NP_/XP_/WP_/YP_/AP_) are
recognized so the identity resolver can label curated vs predicted records.

Parsing is isolated (``parse_fasta``/``parse_esummary``) for offline tests.
"""

from __future__ import annotations

import datetime
import os
import re
import urllib.parse

from .http_util import FetchError, fetch_json, fetch_text

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# RefSeq protein accession prefixes → curation status.
_REFSEQ_CLASS = {
    "NP": "refseq-curated", "XP": "refseq-predicted", "WP": "refseq-nonredundant",
    "YP": "refseq-curated", "AP": "refseq-curated",
}
_REFSEQ_RE = re.compile(r"^(NP|XP|WP|YP|AP)_\d+(\.\d+)?$", re.IGNORECASE)


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _api_key_suffix() -> str:
    key = os.environ.get("NCBI_API_KEY")
    return f"&api_key={urllib.parse.quote(key)}" if key else ""


def refseq_class(accession: str) -> str | None:
    m = _REFSEQ_RE.match((accession or "").strip())
    return _REFSEQ_CLASS.get(m.group(1).upper()) if m else None


def _provenance(acc: str, category: str) -> dict:
    return {
        "source": "NCBI Protein",
        "source_id": acc,
        "source_version": None,
        "retrieved_utc": _now(),
        "evidence_category": category,
        "limitations": ("computationally predicted (XP_)"
                        if refseq_class(acc) == "refseq-predicted" else None),
    }


# ---------------------------------------------------------------------------
# Network entry points
# ---------------------------------------------------------------------------
def fetch_protein(accession: str) -> dict:
    """Fetch a protein record: FASTA sequence + esummary metadata, normalized."""
    acc = (accession or "").strip()
    if not acc:
        raise FetchError("empty NCBI accession")
    quoted = urllib.parse.quote(acc)
    fasta = fetch_text(
        f"{_EUTILS}/efetch.fcgi?db=protein&id={quoted}&rettype=fasta&retmode=text"
        + _api_key_suffix())
    seq_rec = parse_fasta(fasta)
    summary = {}
    try:
        data = fetch_json(
            f"{_EUTILS}/esummary.fcgi?db=protein&id={quoted}&retmode=json"
            + _api_key_suffix())
        summary = parse_esummary(data)
    except FetchError:
        pass  # sequence alone is still useful
    return build_fragment(acc, seq_rec, summary)


# ---------------------------------------------------------------------------
# Pure parsing
# ---------------------------------------------------------------------------
def parse_fasta(text: str) -> dict:
    """Parse a single-record FASTA into {accession, description, organism, sequence}."""
    lines = (text or "").splitlines()
    header = ""
    seq_lines = []
    for line in lines:
        if line.startswith(">"):
            if header:  # stop at a second record
                break
            header = line[1:].strip()
        elif header:
            seq_lines.append(line.strip())
    if not header:
        raise FetchError("no FASTA record in NCBI response")
    parts = header.split(None, 1)
    accession = parts[0] if parts else None
    description = parts[1] if len(parts) > 1 else None
    organism = None
    if description:
        m = re.search(r"\[([^\]]+)\]\s*$", description)
        if m:
            organism = m.group(1)
    return {
        "accession": accession,
        "description": description,
        "organism": organism,
        "sequence": "".join(seq_lines),
    }


def parse_esummary(data: dict) -> dict:
    """Extract the useful fields from an E-utilities esummary JSON payload."""
    result = (data or {}).get("result") or {}
    uids = result.get("uids") or []
    if not uids:
        return {}
    rec = result.get(uids[0]) or {}
    return {
        "title": rec.get("title"),
        "organism": rec.get("organism"),
        "taxid": rec.get("taxid"),
        "length": rec.get("slen"),
        "accession_version": rec.get("accessionversion"),
        "gi": rec.get("gi"),
    }


def build_fragment(acc: str, seq_rec: dict, summary: dict) -> dict:
    cls = refseq_class(acc)
    category = "predicted" if cls == "refseq-predicted" else "curated"
    fragment = {
        "accession": seq_rec.get("accession") or acc,
        "refseq_class": cls,
        "description": seq_rec.get("description"),
        "organism": summary.get("organism") or seq_rec.get("organism"),
        "taxon_id": summary.get("taxid"),
        "sequence": {
            "value": seq_rec.get("sequence"),
            "length": summary.get("length") or len(seq_rec.get("sequence") or ""),
        },
        "title": summary.get("title"),
    }
    return {"data": fragment, "provenance": _provenance(acc, category), "status": "ok"}
