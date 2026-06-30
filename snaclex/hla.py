"""HLA / MHC peptide-binding-groove analysis.

The unifying clinical thread for SnaCleX: from an HLA allele or variant, to the
HLA protein, to a *structural* read-out of how the change affects the
peptide-binding groove. This module:

  1. Detects whether a loaded structure is an HLA/MHC complex (class I or II),
     from chain composition plus optional title/UniProt signals.
  2. Derives the groove as the HLA residues that contact the bound peptide
     (reusing ``interactions.profile_interface`` — the peptide is just a short
     chain), and flags the class-I anchor pockets (B-pocket at peptide P2,
     F-pocket at the C-terminal anchor PΩ).
  3. Classifies mapped variants as anchor-pocket / groove-lining / peripheral,
     a transparent "likely alters the peptide repertoire vs peripheral" call.

Scope: structural interpretation only. This is NOT peptide-binding-affinity or
immunogenicity prediction (use NetMHCpan-style tools for that), and NOT a
clinical or HLA-typing result. Research-only.
"""

from __future__ import annotations

import re

from . import interactions

# Canonical UniProt accessions for HLA heavy chains and class-II subunits.
HLA_UNIPROT = {
    # Class I heavy chains
    "P04439": "HLA-A", "P01889": "HLA-B", "P10321": "HLA-C",
    "P13746": "HLA-A", "P30443": "HLA-A", "P30460": "HLA-B", "P30685": "HLA-B",
    "P30501": "HLA-C", "P03989": "HLA-B", "P10319": "HLA-B",
    "P13747": "HLA-E", "P30511": "HLA-F", "P17693": "HLA-G",
    # beta-2 microglobulin
    "P61769": "B2M",
    # Class II
    "P01903": "HLA-DRA", "P01911": "HLA-DRB1", "P04229": "HLA-DRB1",
    "P01906": "HLA-DQA2", "P01909": "HLA-DQA1", "P01920": "HLA-DQB1",
    "P20036": "HLA-DPA1", "P04440": "HLA-DPB1",
}

_TITLE_RE = re.compile(
    r"\bHLA\b|histocompatibility|\bMHC\b|class\s+I{1,2}\b|class\s+(?:one|two)\b",
    re.IGNORECASE,
)

# Peptide-length windows (residues) for a bound antigen.
_PEPTIDE_MIN, _PEPTIDE_MAX = 7, 25
# Class-I heavy chain and beta-2-microglobulin length windows.
_HEAVY_MIN = 250
_B2M_MIN, _B2M_MAX = 90, 110
# Class-II subunit length window.
_CLASS2_MIN, _CLASS2_MAX = 160, 240


def _chain_residue_counts(structure) -> dict:
    by_chain: dict[str, set] = {}
    for a in structure.protein_atoms:
        by_chain.setdefault(a.chain, set()).add(a.res_seq)
    return {c: len(s) for c, s in by_chain.items()}


def detect(structure, title: str | None = None, uniprots=None) -> dict:
    """Classify a structure as HLA class I / II (or not), with the key chains.

    Detection is geometric (chain lengths) with title/UniProt as corroborating
    signals; either path can establish HLA identity.
    """
    lengths = _chain_residue_counts(structure)
    signals: list[str] = []

    uni_hits = [HLA_UNIPROT[u] for u in (uniprots or []) if u in HLA_UNIPROT]
    if uni_hits:
        signals.append("UniProt: " + ", ".join(sorted(set(uni_hits))))
    if title and _TITLE_RE.search(title):
        signals.append("title mentions HLA/MHC")

    peptides = sorted(
        (c for c, n in lengths.items() if _PEPTIDE_MIN <= n <= _PEPTIDE_MAX),
        key=lambda c: lengths[c],
    )
    heavy = [c for c, n in lengths.items() if n >= _HEAVY_MIN]
    b2m = [c for c, n in lengths.items() if _B2M_MIN <= n <= _B2M_MAX]
    class2 = [c for c, n in lengths.items() if _CLASS2_MIN <= n <= _CLASS2_MAX]

    peptide_chain = peptides[0] if peptides else None
    # HLA gene symbols implied by the mapped UniProt accessions (e.g. HLA-A),
    # excluding beta-2-microglobulin — used for expression-context lookups.
    genes = sorted({HLA_UNIPROT[u] for u in (uniprots or [])
                    if u in HLA_UNIPROT and HLA_UNIPROT[u] != "B2M"})
    result = {
        "is_hla": False,
        "mhc_class": None,
        "confidence": "low",
        "signals": signals,
        "chain_lengths": lengths,
        "heavy_chain": None,
        "b2m_chain": None,
        "class2_chains": [],
        "peptide_chain": peptide_chain,
        "genes": genes,
    }

    # Class I: a long heavy chain (+ optionally beta-2-microglobulin).
    if heavy:
        result["heavy_chain"] = max(heavy, key=lambda c: lengths[c])
        if b2m:
            result["b2m_chain"] = max(b2m, key=lambda c: lengths[c])
        geometric = bool(peptide_chain) and (result["b2m_chain"] is not None)
        if geometric or uni_hits or signals:
            result["is_hla"] = True
            result["mhc_class"] = "I"
            result["confidence"] = "high" if (geometric and signals) else (
                "medium" if geometric or uni_hits else "low"
            )
            return result

    # Class II: two medium subunits (+ a bound peptide).
    if len(class2) >= 2:
        result["class2_chains"] = sorted(class2, key=lambda c: -lengths[c])[:2]
        geometric = bool(peptide_chain)
        if geometric or uni_hits or signals:
            result["is_hla"] = True
            result["mhc_class"] = "II"
            result["confidence"] = "high" if (geometric and signals) else "medium"
            return result

    # Signal-only fallback (e.g. apo HLA with no peptide modeled).
    if uni_hits:
        result["is_hla"] = True
        result["confidence"] = "low"
    return result


def _hla_chains(detection: dict) -> list:
    if detection.get("mhc_class") == "II":
        return list(detection.get("class2_chains") or [])
    heavy = detection.get("heavy_chain")
    return [heavy] if heavy else []


def analyze_groove(structure, detection: dict) -> dict:
    """Derive groove + anchor-pocket residues from peptide contacts."""
    peptide = detection.get("peptide_chain")
    hla_chains = _hla_chains(detection)
    if not peptide or not hla_chains:
        return {
            "available": False,
            "reason": "No bound peptide and HLA chain pair were identified "
                      "(an apo structure, or peptide not modeled).",
        }

    prof = interactions.profile_interface(structure, hla_chains, [peptide])
    groove_residues = prof["interface_residues_a"]
    peptide_residues = prof["interface_residues_b"]

    # Full peptide numbering (N->C), so anchors are true P2 / C-terminus even if
    # some positions make no HLA contact. Class-I anchors are P2 (B-pocket) and
    # the C-terminal residue PΩ (F-pocket).
    all_pep_seqs = sorted(
        {a.res_seq for a in structure.protein_atoms if a.chain == peptide}
    )
    anchor_pep_seqs: set[int] = set()
    if detection.get("mhc_class") == "I" and len(all_pep_seqs) >= 2:
        anchor_pep_seqs = {all_pep_seqs[1], all_pep_seqs[-1]}

    # HLA residues contacting an anchor peptide position = anchor-pocket residues.
    anchor_ids: set[str] = set()
    for ix in prof["interactions"]:
        if ix["atom_b"]["res_seq"] in anchor_pep_seqs:
            anchor_ids.add(ix["atom_a"]["res_id"])

    for r in groove_residues:
        r["pocket"] = "anchor" if r["res_id"] in anchor_ids else "groove"

    return {
        "available": True,
        "mhc_class": detection.get("mhc_class"),
        "hla_chains": hla_chains,
        "peptide_chain": peptide,
        "peptide_length": len(all_pep_seqs),
        "groove_residues": groove_residues,
        "groove_residue_ids": [r["res_id"] for r in groove_residues],
        "anchor_pocket_residue_ids": sorted(anchor_ids),
        "peptide_residues": peptide_residues,
        "anchor_peptide_positions": sorted(anchor_pep_seqs),
        "counts": prof["counts"],
        "interaction_total": prof["interaction_total"],
    }


def classify_variants(groove: dict, mapped_variants) -> list:
    """Tag each mapped variant by its location in the groove.

    anchor-pocket  -> directly shapes an anchor pocket; most likely to alter the
                      bound-peptide repertoire.
    groove-lining  -> lines the cleft; may modulate peptide presentation.
    peripheral     -> outside the modeled peptide contacts.
    """
    if not groove.get("available"):
        return []
    groove_ids = set(groove.get("groove_residue_ids") or [])
    anchor_ids = set(groove.get("anchor_pocket_residue_ids") or [])
    out = []
    for v in mapped_variants:
        if not v.get("mapped"):
            continue
        rid = v.get("res_id")
        if rid in anchor_ids:
            cls, note = "anchor-pocket", (
                "Shapes an anchor pocket — most likely to alter which peptides "
                "the allele can present."
            )
        elif rid in groove_ids:
            cls, note = "groove-lining", (
                "Lines the peptide-binding groove — may modulate presentation."
            )
        else:
            cls, note = "peripheral", (
                "Outside the modeled peptide contacts in this structure."
            )
        out.append({
            "input": v.get("input"),
            "res_id": rid,
            "wt": v.get("wt"),
            "mut": v.get("mut"),
            "groove_class": cls,
            "interpretation": note,
        })
    return out


def curated_cases() -> list:
    """Well-known HLA structures for examples / benchmarking (fetched live)."""
    return [
        {"allele": "HLA-A*02:01", "pdb": "1HHK", "mhc_class": "I",
         "note": "Common class-I allele; central to neoantigen presentation."},
        {"allele": "HLA-B*57:01", "pdb": "3VRI", "mhc_class": "I",
         "note": "Abacavir hypersensitivity — drug alters the bound peptide repertoire."},
        {"allele": "HLA-B*27:05", "pdb": "1HSA", "mhc_class": "I",
         "note": "Associated with ankylosing spondylitis (autoimmunity)."},
        {"allele": "HLA-DQ2.5", "pdb": "1S9V", "mhc_class": "II",
         "note": "Presents deamidated gliadin peptides in celiac disease."},
    ]
