"""Ensembl VEP: genomic variant -> protein consequence (optional, env-gated).

Turns a genomic SNV (``chr:pos ref>alt``) into its protein consequence — gene,
UniProt accession, protein position, and amino-acid change — via the Ensembl
VEP REST API. That lets the server convert genomic input into a protein-position
variant and map it onto the structure through the existing residue mapper.

OFF unless ``SNACLEX_ENABLE_VEP`` is set, since it adds an external dependency
and the remote env's network policy may block Ensembl. The HTTP call sits behind
an injectable seam so orchestration is testable offline, and any failure degrades
to ``None`` (the variant simply keeps its TOPMed/BRAVO frequency without a
structural location). Research-only: a consequence annotation, not a clinical call.
"""

from __future__ import annotations

import os

from .http_util import FetchError, fetch_json

_ENV_FLAG = "SNACLEX_ENABLE_VEP"
_BASE = "https://rest.ensembl.org"


def available() -> bool:
    """True only when VEP annotation is explicitly enabled."""
    v = (os.environ.get(_ENV_FLAG) or "").strip().lower()
    return v not in ("", "0", "off", "false", "no")


def _region_url(chrom, pos, ref, alt) -> str:
    chrom = str(chrom).replace("chr", "")
    end = int(pos) + len(str(ref)) - 1
    region = f"{chrom}:{pos}-{end}"
    return (f"{_BASE}/vep/human/region/{region}/{alt}"
            "?content-type=application/json&uniprot=1&canonical=1")


def _strip_version(acc):
    """'P04637.4' -> 'P04637'; passthrough for already-bare accessions."""
    if not acc:
        return None
    if isinstance(acc, list):
        acc = acc[0] if acc else None
    return str(acc).split(".")[0] if acc else None


def _pick_consequence(payload) -> dict | None:
    """Pick the best protein-coding consequence (canonical + has a UniProt acc).

    Returns ``{gene, uniprot, protein_position, wt_aa, mut_aa, consequence}`` or
    None when no protein-altering, UniProt-mapped consequence is present.
    """
    if isinstance(payload, list):
        record = payload[0] if payload else None
    else:
        record = payload
    if not isinstance(record, dict):
        return None
    cons = record.get("transcript_consequences") or []
    # Prefer canonical transcripts, then those carrying a UniProt (swissprot) id.
    cons = sorted(
        cons,
        key=lambda c: (0 if c.get("canonical") else 1, 0 if c.get("swissprot") else 1),
    )
    for c in cons:
        aa = c.get("amino_acids")          # e.g. "R/H"
        pos = c.get("protein_start")
        uni = _strip_version(c.get("swissprot"))
        if not aa or pos is None or not uni or "/" not in aa:
            continue
        wt, mut = aa.split("/", 1)
        if not wt or not mut or wt == mut:
            continue
        return {
            "gene": c.get("gene_symbol"),
            "uniprot": uni,
            "protein_position": int(pos),
            "wt_aa": wt,
            "mut_aa": mut,
            "consequence": (c.get("consequence_terms") or [None])[0],
        }
    return None


def annotate_one(chrom, pos, ref, alt, fetcher=None) -> dict | None:
    """Return the protein consequence for one SNV, or None (never raises)."""
    fetcher = fetcher or fetch_json
    try:
        data = fetcher(_region_url(chrom, pos, ref, alt))
    except FetchError:
        return None
    return _pick_consequence(data)
