"""TOPMed/BRAVO public allele-frequency client (supplemental, best-effort).

BRAVO (bravo.sph.umich.edu) serves *aggregate* allele frequencies computed from
TOPMed. Only public aggregate frequencies are used here — never individual
genotypes or controlled-access data (those stay in dbGaP / BioData Catalyst /
Terra). Many BRAVO deployments gate queries behind login, so this client is
best-effort: it returns a structured ``{"available": False, ...}`` rather than
raising, and the base URL is overridable via ``SNACLEX_BRAVO_API`` for sites
that expose an open endpoint or a mirror.

Frequencies are keyed by genomic coordinate, so this is reached only on the
genomic (VEP) variant path, not on protein-position input.
"""

from __future__ import annotations

import os

from .http_util import FetchError, fetch_json

_DEFAULT_BASE = "https://bravo.sph.umich.edu/freeze8/hg38/api"


def _base() -> str:
    return (os.environ.get("SNACLEX_BRAVO_API") or _DEFAULT_BASE).rstrip("/")


def _unavailable(reason: str) -> dict:
    return {"available": False, "reason": reason}


def _coerce_freq(data) -> dict | None:
    """Pull (allele_freq, allele_count, allele_num) out of a BRAVO-style payload.

    BRAVO responses vary by deployment/freeze; we look for the common keys and
    fall back gracefully. Returns None if no frequency field is present.
    """
    if isinstance(data, dict) and "data" in data and isinstance(data["data"], list):
        records = data["data"]
        rec = records[0] if records else None
    elif isinstance(data, list):
        rec = data[0] if data else None
    else:
        rec = data
    if not isinstance(rec, dict):
        return None

    af = rec.get("allele_freq")
    if af is None:
        af = rec.get("af")
    ac = rec.get("allele_count", rec.get("ac"))
    an = rec.get("allele_num", rec.get("an"))
    if af is None and ac is not None and an:
        try:
            af = ac / an
        except (TypeError, ZeroDivisionError):
            af = None
    if af is None:
        return None
    return {
        "allele_freq": af,
        "allele_count": ac,
        "allele_num": an,
        "rsid": rec.get("rsids") or rec.get("rsid"),
    }


def variant_frequency(chrom, pos, ref, alt, base: str | None = None) -> dict:
    """Look up the TOPMed allele frequency for one SNV (best-effort).

    Returns ``{"available": True, "allele_freq": ..., ...}`` on success, or
    ``{"available": False, "reason": ...}`` if the lookup is unavailable.
    """
    chrom = str(chrom).replace("chr", "")
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    url = f"{base or _base()}/variant/api/snv/{variant_id}"
    try:
        data = fetch_json(url)
    except FetchError as exc:
        return _unavailable(f"BRAVO lookup unavailable ({exc})")
    freq = _coerce_freq(data)
    if freq is None:
        return _unavailable("variant not found in TOPMed/BRAVO (or no frequency field)")
    out = {"available": True, "variant_id": variant_id, "source": "TOPMed/BRAVO"}
    out.update(freq)
    return out


def classify_frequency(af) -> str:
    """Bucket an allele frequency into a coarse, human-readable rarity label."""
    if af is None:
        return "unknown"
    if af == 0:
        return "absent"
    if af < 0.0001:
        return "ultra-rare (<0.01%)"
    if af < 0.01:
        return "rare (<1%)"
    if af < 0.05:
        return "low-frequency (1-5%)"
    return "common (>5%)"
