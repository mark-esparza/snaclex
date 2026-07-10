"""Variant / genetics analysis (Phase 3 — Genetics workspace).

Turns a variant input (``BRAF V600E``, ``TP53 R175H``, ``P15056:p.Val600Glu``)
into a residue-level analysis against a ProteinRecord: sequence-coordinate
mapping with a **WT-residue validation** (catches wrong-isoform / off-by-one
numbering), coding consequence, domain-disruption context, overlap with curated
residue annotations and known natural variants, and any clinical interpretation.

Scientific safeguard (non-negotiable): a variant is **never** interpreted
clinically. The analysis displays the originating evidence, the record's review
status, and explicit limitations — it makes no diagnosis, prognosis, or
recommendation. Parsing and mapping are pure and offline-testable.
"""

from __future__ import annotations

import re

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "SEC": "U", "PYL": "O", "TER": "*", "STOP": "*",
}
_AA1 = set("ACDEFGHIKLMNPQRSTVWYUO*")

_HGVS_RE = re.compile(
    r"^(?P<acc>\S+):p\.(?P<wt>[A-Za-z]{1,3})(?P<pos>\d+)(?P<mut>[A-Za-z]{1,3}|\*|=)$")
_SIMPLE_RE = re.compile(
    r"^(?P<gene>\S+)\s+(?P<wt>[A-Za-z]{1,3})(?P<pos>\d+)(?P<mut>[A-Za-z]{1,3}|\*)$")

_DISCLAIMER = (
    "Not a clinical interpretation. This displays the originating database "
    "evidence, its review status, and limitations only; it makes no diagnosis, "
    "prognosis, or treatment recommendation. Confirm with an appropriate "
    "clinical resource and review status."
)


class VariantError(ValueError):
    """Raised when a variant string cannot be parsed."""


def _norm_aa(token: str) -> str:
    t = token.strip().upper()
    if t in ("=",):
        return "="  # HGVS synonymous marker
    if len(t) == 1:
        if t in _AA1:
            return t
        raise VariantError(f"unknown amino-acid code {token!r}")
    if t in _AA3TO1:
        return _AA3TO1[t]
    raise VariantError(f"unknown amino-acid code {token!r}")


def parse(query: str) -> dict:
    """Parse a variant string into {gene_or_acc, wt, position, mut, notation}."""
    q = (query or "").strip()
    m = _HGVS_RE.match(q)
    kind = "hgvs"
    if not m:
        m = _SIMPLE_RE.match(q)
        kind = "simple"
    if not m:
        raise VariantError(
            f"could not parse variant {query!r} (expected e.g. 'BRAF V600E' or "
            "'P15056:p.Val600Glu')")
    gene_or_acc = m.group("acc") if kind == "hgvs" else m.group("gene")
    wt = _norm_aa(m.group("wt"))
    mut = "=" if m.group("mut") == "=" else _norm_aa(m.group("mut"))
    if mut == "=":
        mut = wt
    pos = int(m.group("pos"))
    return {
        "gene_or_acc": gene_or_acc, "wt": wt, "position": pos, "mut": mut,
        "notation": f"{wt}{pos}{mut}", "input": q,
    }


def consequence(wt: str, mut: str) -> str:
    """Coarse coding consequence at the protein level.

    Frameshift / splice / start-loss / stop-loss need transcript (nucleotide)
    context and are out of scope here (they require Ensembl/VEP, Phase 3).
    """
    if mut == "*":
        return "nonsense (stop gain)"
    if wt == mut:
        return "synonymous"
    return "missense"


def _covers(feature: dict, pos: int) -> bool:
    start = feature.get("start")
    end = feature.get("end", start)
    if start is None:
        return False
    if end is None:
        end = start
    return start <= pos <= end


def analyze(record: dict, variant: dict) -> dict:
    """Map a parsed variant onto a ProteinRecord (pure, offline-testable)."""
    seq = (record.get("sequence") or {}).get("value") or ""
    length = len(seq)
    pos = variant["position"]
    wt_expected = variant["wt"]

    warnings: list[str] = []
    wt_in_sequence = seq[pos - 1] if 1 <= pos <= length else None
    in_range = wt_in_sequence is not None
    wt_matches = in_range and wt_in_sequence == wt_expected
    if not in_range:
        warnings.append(
            f"position {pos} is outside the canonical sequence (length {length}) — "
            "check isoform/numbering")
    elif not wt_matches:
        warnings.append(
            f"WT residue mismatch: variant says {wt_expected}{pos} but the canonical "
            f"sequence has {wt_in_sequence}{pos} — likely a different isoform or "
            "numbering scheme")

    domains_affected = [
        {"name": d.get("name"), "type": d.get("type"),
         "start": d.get("start"), "end": d.get("end")}
        for d in (record.get("domains_motifs") or []) if _covers(d, pos)]

    residue_context = [
        {"kind": r.get("kind"), "description": r.get("description")}
        for r in (record.get("residue_annotations") or []) if _covers(r, pos)]
    ptm_context = [
        {"type": p.get("type"), "description": p.get("description")}
        for p in (record.get("ptms") or []) if p.get("position") == pos]

    known = next((v for v in (record.get("variants") or [])
                  if v.get("position") == pos), None)

    review_status = record.get("review_status")
    interpretations = []
    if known and known.get("description"):
        interpretations.append({"source": "UniProt (natural variant)",
                                "text": known["description"]})
    for dis in record.get("disease_associations") or []:
        if dis.get("name"):
            interpretations.append({"source": "UniProt (disease)",
                                    "text": dis.get("description") or dis["name"],
                                    "name": dis.get("name")})

    structures = record.get("structures") or {}
    experimental = [s.get("pdb_id") for s in structures.get("experimental", [])]

    return {
        "notation": variant["notation"],
        "consequence": consequence(wt_expected, variant["mut"]),
        "sequence_mapping": {
            "canonical_position": pos,
            "wt_expected": wt_expected,
            "wt_in_sequence": wt_in_sequence,
            "wt_matches": wt_matches,
            "sequence_length": length,
        },
        "domains_affected": domains_affected,
        "domain_disruption": bool(domains_affected),
        "residue_context": residue_context + ptm_context,
        "known_variant": known,
        "clinical": {
            "interpretations": interpretations,
            "review_status": review_status,
            "disclaimer": _DISCLAIMER,
        },
        "structure_mapping": {
            "experimental_structures": experimental,
            "note": ("Canonical (UniProt) coordinate given; per-PDB author residue "
                     "numbering requires the SIFTS / RCSB Sequence Coordinates "
                     "mapping (Structural Analysis)."),
        },
        "limitations": [
            "Protein-level consequence only; frameshift/splice/start-loss/stop-loss "
            "need transcript (nucleotide) context (Ensembl/VEP).",
            "Conservation at this residue is available for structures via "
            "/api/evolution; not computed here.",
            "Population frequency and full clinical review require an authorized "
            "source (e.g. ClinVar) — not integrated in this slice.",
        ],
        "warnings": warnings,
        "method": {"tool": "SnaCleX variant mapping",
                   "note": "Sequence-coordinate mapping + curated-annotation overlap; "
                           "no clinical interpretation is performed."},
    }
