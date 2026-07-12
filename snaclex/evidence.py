"""Universal evidence model + A–F evidence levels (Stages 6 & 7).

Every "protein X does Y to chemical/protein Z" claim becomes one evidence object
at exactly one level; levels are never flattened into an undifferentiated
"interaction". The module enforces the hard scientific invariants from
docs/platform/06-evidence-ranking.md:

* A docking score is a *fit score* (``units == "score"``), never an affinity.
* Homology-transferred evidence (Level D) must carry its transfer basis.
* ML/similarity evidence (Level F) must carry its computation/transfer basis.
* A pair's evidence is summarized as per-level *counts*, never a boolean.

Builders (``level_a_from_pdb_ligand`` etc.) are pure and offline-testable.
"""

from __future__ import annotations

import datetime
import hashlib

# Assay/units that must never appear on a docking (Level E) object.
_AFFINITY_ASSAYS = {"KI", "KD", "IC50", "EC50", "AC50", "KM"}
_AFFINITY_UNITS = {"NM", "UM", "µM", "PM", "M", "MM", "NMOL"}

LEVELS = {
    "A": {"name": "Direct experimental structural", "category": "experimental", "direct": True},
    "B": {"name": "Direct biochemical/cellular", "category": "biochemical", "direct": True},
    "C": {"name": "Curated database association", "category": "curated", "direct": False},
    "D": {"name": "Homology-transferred", "category": "homology", "direct": False},
    "E": {"name": "Computational docking", "category": "docking", "direct": False},
    "F": {"name": "ML / similarity prediction", "category": "ml", "direct": False},
}


class EvidenceInvariantError(ValueError):
    """Raised when an evidence object violates a scientific invariant."""


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _evidence_id(subject_ref, predicate, object_ref, source_record) -> str:
    raw = f"{subject_ref}|{predicate}|{object_ref}|{source_record}"
    return "SNX:EV:" + hashlib.sha256(raw.encode()).hexdigest()[:16]


def make_evidence(*, subject_ref, predicate, object_ref, level, source,
                  source_record=None, assay_type=None, value=None, units=None,
                  relation=None, outcome=None, experimental_method=None,
                  species=None, protein_construct=None, sequence_or_isoform=None,
                  chemical_form=None, publication=None, confidence=None,
                  curation_status=None, computation=None, transfer=None,
                  limitations=None) -> dict:
    """Construct and validate a universal evidence object."""
    if level not in LEVELS:
        raise EvidenceInvariantError(f"unknown evidence level {level!r}")
    meta = LEVELS[level]

    # --- invariants -------------------------------------------------------
    if level == "E":
        if (assay_type or "").upper() in _AFFINITY_ASSAYS:
            raise EvidenceInvariantError("docking (E) may not carry an affinity assay type")
        if units and units.upper() in _AFFINITY_UNITS:
            raise EvidenceInvariantError("docking (E) may not carry affinity units")
        units = "score"
        if not computation:
            raise EvidenceInvariantError("docking (E) must carry a computation block")
    if level == "D" and not (transfer and transfer.get("sequence_identity") is not None):
        raise EvidenceInvariantError("homology transfer (D) must carry transfer.sequence_identity")
    if level == "F" and not (computation or transfer):
        raise EvidenceInvariantError("ML/similarity (F) must carry a computation or transfer basis")

    return {
        "evidence_id": _evidence_id(subject_ref, predicate, object_ref, source_record),
        "subject": subject_ref,
        "predicate": predicate,
        "object": object_ref,
        "level": level,
        "level_name": meta["name"],
        "category": meta["category"],
        "is_direct": meta["direct"],
        "source": source,
        "source_record": source_record,
        "experimental_method": experimental_method,
        "assay_type": assay_type,
        "value": value,
        "units": units,
        "relation": relation,
        "outcome": outcome,
        "species": species,
        "protein_construct": protein_construct,
        "sequence_or_isoform": sequence_or_isoform,
        "chemical_form": chemical_form,
        "publication": publication,
        "confidence": confidence,
        "curation_status": curation_status,
        "computation": computation,
        "transfer": transfer,
        "limitations": limitations,
        "retrieved_utc": _now(),
    }


# ---------------------------------------------------------------------------
# Stage 7 — chemical-form match classification
# ---------------------------------------------------------------------------
def classify_chemical_match(query: dict, candidate: dict) -> str:
    """Classify how a candidate chemical relates to the queried chemical.

    Returns one of: exact | parent | stereo | salt | substructure | name | unresolved.
    ``query``/``candidate`` are dicts with optional cid, inchikey, parent_cid.
    Salts/parents/stereoisomers are never silently treated as identical.
    """
    qk = (query.get("inchikey") or "")
    ck = (candidate.get("inchikey") or "")
    if query.get("cid") and query.get("cid") == candidate.get("cid"):
        return "exact"
    if qk and ck:
        if qk == ck:
            return "exact"
        # InChIKey: block1 = skeleton, block2 = stereo/isotope/protonation layer.
        if qk[:14] == ck[:14]:
            return "stereo"  # same connectivity skeleton, differing stereo/protonation
    if query.get("parent_cid") and candidate.get("cid") == query.get("parent_cid"):
        return "parent"
    if candidate.get("parent_cid") and candidate.get("parent_cid") == query.get("cid"):
        return "salt"
    if candidate.get("name_only"):
        return "name"
    return "unresolved"


# ---------------------------------------------------------------------------
# Evidence builders
# ---------------------------------------------------------------------------
def level_a_from_pdb_ligand(*, protein_ref, chemical_ref, pdb_id, chain,
                            ligand_id, chemical_form=None, contacts=None,
                            method=None, resolution_a=None,
                            protein_construct=None, species=None) -> dict:
    """Level A: a chemical bound to the requested protein in a PDB structure."""
    return make_evidence(
        subject_ref=protein_ref, predicate="binds", object_ref=chemical_ref,
        level="A", source="RCSB PDB",
        source_record=f"PDB {pdb_id} chain {chain} ligand {ligand_id}",
        experimental_method=method, species=species,
        protein_construct=protein_construct, chemical_form=chemical_form,
        confidence=None, curation_status="experimental",
        computation={"contacts": contacts, "resolution_A": resolution_a} if contacts else None,
        limitations="crystallographic construct may differ from the canonical sequence")


def level_b_from_bioassay(*, protein_ref, chemical_ref, assay_row,
                          chemical_form=None) -> dict:
    """Level B: a PubChem BioAssay row measuring the compound (incl. inactives)."""
    return make_evidence(
        subject_ref=protein_ref, predicate="inhibits", object_ref=chemical_ref,
        level="B", source="PubChem BioAssay",
        source_record=f"AID {assay_row.get('aid')}",
        assay_type=assay_row.get("activity_name"),
        value=assay_row.get("activity_value_uM"),
        units="uM" if assay_row.get("activity_value_uM") is not None else None,
        outcome=assay_row.get("outcome"),
        experimental_method=assay_row.get("assay_type"),
        publication=(f"PMID:{assay_row['pmid']}" if assay_row.get("pmid") else None),
        chemical_form=chemical_form, curation_status="curated",
        limitations="assay target definition may not match this exact isoform/construct")


def level_d_transfer(*, protein_ref, chemical_ref, from_protein,
                     sequence_identity, aligned_site_identity=None,
                     species=None, source_record=None) -> dict:
    """Level D: compound binds a homolog; evidence transferred with its basis."""
    return make_evidence(
        subject_ref=protein_ref, predicate="binds", object_ref=chemical_ref,
        level="D", source="homology transfer", source_record=source_record,
        species=species,
        transfer={"from_protein": from_protein, "sequence_identity": sequence_identity,
                  "aligned_site_identity": aligned_site_identity,
                  "limitations": "binding inferred from a homolog, not observed on this protein"},
        confidence=sequence_identity,
        limitations="inference from a related protein — not direct evidence for this protein")


def docking_evidence(*, protein_ref, chemical_ref, score, software, params,
                     seed=None, structure_ref=None, per_atom_score=None) -> dict:
    """Level E: a docking pose fit score (never an affinity)."""
    return make_evidence(
        subject_ref=protein_ref, predicate="binds", object_ref=chemical_ref,
        level="E", source=software, source_record=structure_ref,
        value=score, units="score", relation="=",
        computation={"software": software, "params": params,
                     "reproducibility": {"seed": seed},
                     "per_atom_score": per_atom_score},
        limitations=("predicted fit score (lower = better geometric fit); NOT a "
                     "binding affinity, NOT kcal/mol, NOT experimental validation"))


# ---------------------------------------------------------------------------
# Summaries (never a boolean)
# ---------------------------------------------------------------------------
def summarize(evidence_list: list[dict]) -> dict:
    """Per-level counts + direct-vs-predicted split (no 'interacts' boolean)."""
    counts = {lvl: 0 for lvl in LEVELS}
    for e in evidence_list:
        lvl = e.get("level")
        if lvl in counts:
            counts[lvl] += 1
    direct = counts["A"] + counts["B"]
    predicted = counts["C"] + counts["D"] + counts["E"] + counts["F"]
    return {"by_level": counts, "direct_experimental": direct,
            "predicted_or_inferred": predicted,
            "note": "Counts by evidence level. Direct (A/B) and predicted/inferred "
                    "(C–F) are reported separately and never merged into a single claim."}


def rank(evidence_list: list[dict]) -> list[dict]:
    """Order evidence by directness A>B>C>D>E>F (stable within a level)."""
    order = {lvl: i for i, lvl in enumerate("ABCDEF")}
    return sorted(evidence_list, key=lambda e: order.get(e.get("level"), 99))
