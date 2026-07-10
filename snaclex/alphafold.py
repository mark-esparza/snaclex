"""AlphaFold DB adapter — the predicted-structure layer.

Retrieves database-predicted models for a UniProt accession and exposes model id,
UniProt mapping, coverage, fragmentation, version/source, per-residue pLDDT
confidence, and PAE availability. An AlphaFold model is **never** presented as
experimental — the record builder badges it "predicted" and the docking gate
(Stage 4/6) refuses to auto-dock low-confidence regions.

pLDDT is stored in the B-factor column of the model coordinates; ``parse_plddt``
reads it directly from PDB ATOM records so confidence bands and low-confidence
regions are derived from the real model, not guessed. ``globalMetricValue`` from
the API is treated as the mean-pLDDT summary (marked VERIFY until round-tripped).
"""

from __future__ import annotations

import datetime

from .http_util import FetchError, fetch_json, fetch_text

_API = "https://alphafold.ebi.ac.uk/api/prediction"

# pLDDT confidence bands (AlphaFold's documented thresholds).
_BANDS = (("very_high", 90.0), ("confident", 70.0), ("low", 50.0), ("very_low", 0.0))


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _provenance(model_id, version) -> dict:
    return {
        "source": "AlphaFold DB",
        "source_id": model_id,
        "source_version": version,
        "retrieved_utc": _now(),
        "evidence_category": "predicted",
        "limitations": "computationally predicted model — not an experimental structure",
    }


# ---------------------------------------------------------------------------
# Network entry points
# ---------------------------------------------------------------------------
def fetch_models(uniprot_accession: str) -> list[dict]:
    """Return normalized predicted-model records for a UniProt accession.

    Empty list (not an error) when AlphaFold DB has no model — the caller then
    falls through to sequence-only analysis.
    """
    acc = (uniprot_accession or "").strip()
    if not acc:
        raise FetchError("empty UniProt accession")
    try:
        data = fetch_json(f"{_API}/{acc}")
    except FetchError:
        return []
    return parse_prediction(data)


def fetch_model_pdb(pdb_url: str) -> str:
    """Download a model's PDB text (for pLDDT extraction / viewing)."""
    return fetch_text(pdb_url)


# ---------------------------------------------------------------------------
# Pure parsing
# ---------------------------------------------------------------------------
def parse_prediction(data) -> list[dict]:
    """Normalize the AlphaFold prediction API list into model records."""
    if not isinstance(data, list):
        return []
    models = []
    for m in data:
        start = m.get("uniprotStart")
        end = m.get("uniprotEnd")
        seq_len = len(m.get("uniprotSequence") or "") or None
        coverage = None
        if start is not None and end is not None and seq_len:
            coverage = round((end - start + 1) / seq_len, 4)
        model_id = m.get("entryId")
        version = m.get("latestVersion")
        models.append({
            "data": {
                "model_id": model_id,
                "uniprot_accession": m.get("uniprotAccession"),
                "uniprot_start": start,
                "uniprot_end": end,
                "coverage_fraction": coverage,
                "fragmented": bool(len(data) > 1),
                "version": version,
                "all_versions": m.get("allVersions"),
                "created": m.get("modelCreatedDate"),
                "global_plddt_mean": m.get("globalMetricValue"),  # VERIFY semantics
                "pdb_url": m.get("pdbUrl"),
                "cif_url": m.get("cifUrl"),
                "bcif_url": m.get("bcifUrl"),
                "pae_available": bool(m.get("paeImageUrl") or m.get("paeDocUrl")),
                "pae_image_url": m.get("paeImageUrl"),
                "pae_doc_url": m.get("paeDocUrl"),
                "kind": "predicted",
                "source": "AlphaFold DB",
            },
            "provenance": _provenance(model_id, version),
            "status": "ok",
        })
    return models


def parse_plddt(pdb_text: str) -> list[float]:
    """Extract per-residue pLDDT (CA B-factor) from AlphaFold PDB text.

    Columns 61-66 hold the B-factor in fixed-width PDB ATOM records; AlphaFold
    stores pLDDT there. Returns one value per residue (CA atom), in order.
    """
    plddt: list[float] = []
    for line in (pdb_text or "").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                plddt.append(float(line[60:66]))
            except ValueError:
                continue
    return plddt


def confidence_bands(plddt: list[float]) -> dict:
    """Fraction of residues in each AlphaFold confidence band."""
    n = len(plddt)
    if not n:
        return {name: None for name, _ in _BANDS}
    counts = {name: 0 for name, _ in _BANDS}
    for v in plddt:
        for name, lo in _BANDS:
            if v >= lo:
                counts[name] += 1
                break
    return {name: round(counts[name] / n, 4) for name, _ in _BANDS}


def low_confidence_regions(plddt: list[float], threshold: float = 70.0,
                           min_len: int = 3) -> list[dict]:
    """Contiguous residue ranges (1-based) below a pLDDT threshold."""
    regions: list[dict] = []
    start = None
    for i, v in enumerate(plddt, start=1):
        if v < threshold:
            if start is None:
                start = i
        else:
            if start is not None and i - start >= min_len:
                regions.append({"start": start, "end": i - 1})
            start = None
    if start is not None and len(plddt) - start + 1 >= min_len:
        regions.append({"start": start, "end": len(plddt)})
    return regions


def docking_assessment(model: dict, plddt: list[float] | None = None) -> dict:
    """Decide whether a predicted model is suitable for docking (Stage 4 gate).

    Never auto-dock a low-confidence model: require mean pLDDT ≥ 70 as a coarse
    gate here; production adds per-pocket confidence + preparation checks.
    """
    mean = None
    if plddt:
        mean = sum(plddt) / len(plddt)
    elif model.get("global_plddt_mean") is not None:
        mean = model["global_plddt_mean"]
    suitable = mean is not None and mean >= 70.0
    reason = None
    if not suitable:
        reason = ("insufficient confidence (mean pLDDT < 70 or unknown); "
                  "requires per-region confidence + preparation assessment")
    return {"docking_suitable": suitable, "mean_plddt": mean,
            "docking_block_reason": reason}
