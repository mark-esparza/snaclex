"""Unified, sequence-first Protein Record (the platform's central aggregate).

A ProteinRecord exists independently of any structure. ``build_record`` composes
already-fetched adapter fragments (so it is offline-testable) into the shape in
docs/platform/04-data-model.md, runs the pure sequence analysis, and stamps
provenance. ``structure_availability`` performs the Stage-4 hierarchy
(experimental PDB → predicted AlphaFold → sequence-only) with a docking gate.

Identity is never inferred from name/gene alone — the canonical id is built from
taxon + accession + isoform role, matching the resolver's identity key.
"""

from __future__ import annotations

import datetime

from . import __version__, alphafold, checksum, homology, rcsb, seqanalysis
from .http_util import FetchError


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def canonical_id(taxon_id, accession, isoform_ordinal=1) -> str:
    return f"SNX:PRT:{taxon_id or 'NA'}:{accession or 'SEQONLY'}:{isoform_ordinal}"


def _sequence_from(resolution: dict, up: dict | None, ncbi_rec: dict | None) -> str | None:
    if up and (up.get("sequence") or {}).get("value"):
        return up["sequence"]["value"]
    if resolution.get("sequence"):
        return resolution["sequence"]
    if ncbi_rec and (ncbi_rec.get("sequence") or {}).get("value"):
        return ncbi_rec["sequence"]["value"]
    return None


def build_record(resolution: dict, *, ph: float = 7.0) -> dict:
    """Assemble a ProteinRecord from a resolver result. Pure (no network)."""
    up = resolution.get("_uniprot_fragment")
    upi = resolution.get("_uniparc_fragment")
    ncbi_rec = resolution.get("ncbi")
    seq = _sequence_from(resolution, up, ncbi_rec)
    warnings = list(resolution.get("warnings") or [])
    sources_used: list[str] = []

    org = (up or {}).get("organism") or {}
    taxon_id = org.get("taxon_id") or (ncbi_rec or {}).get("taxon_id")
    acc = (up or {}).get("accession")
    review = (up or {}).get("review_status") or (
        (ncbi_rec or {}).get("refseq_class") and "ncbi") or "sequence-only"

    if up:
        sources_used.append("UniProtKB")
    if upi:
        sources_used.append("UniParc")
    if ncbi_rec:
        sources_used.append("NCBI Protein")

    # Sequence block + checksums (compute if UniProt did not supply them).
    seq_block = {}
    if up:
        seq_block = dict(up.get("sequence") or {})
    if seq:
        computed = checksum.checksums(seq)
        seq_block.setdefault("value", seq)
        seq_block.setdefault("length", len(seq))
        seq_block["checksum"] = {
            "crc64": seq_block.get("crc64") or computed["crc64"],
            "md5": seq_block.get("md5") or computed["md5"],
            "sha256": computed["sha256"],
        }

    if review == "unreviewed":
        warnings.append("unreviewed record (TrEMBL)")
    if not up and seq:
        warnings.append("sequence-only record — no UniProt annotation resolved")

    # Sequence-only analysis (always available; labelled calculated).
    seq_analysis = None
    if seq:
        try:
            seq_analysis = seqanalysis.analyze(seq, ph=ph)
        except seqanalysis.SequenceError as exc:
            warnings.append(f"sequence analysis skipped: {exc}")

    record = {
        "canonical_id": canonical_id(taxon_id, acc),
        "preferred_name": (up or {}).get("preferred_name")
                          or (ncbi_rec or {}).get("description"),
        "gene": (up or {}).get("gene") or {"symbol": None, "synonyms": []},
        "organism": {"scientific_name": org.get("scientific_name")
                     or (ncbi_rec or {}).get("organism"), "taxon_id": taxon_id},
        "accessions": {
            "uniprot_primary": acc,
            "uniprot_secondary": (up or {}).get("secondary_accessions") or [],
            "uniparc": (upi or {}).get("upi"),
            "ncbi_protein": [ncbi_rec["accession"]] if ncbi_rec and ncbi_rec.get("accession") else [],
            "refseq": [e["id"] for e in ((up or {}).get("cross_references") or {}).get("RefSeq", [])],
            "pdb_entities": [{"pdb": e["id"]} for e in ((up or {}).get("cross_references") or {}).get("PDB", [])],
            "alphafold": [f"AF-{acc}-F1"] if acc else [],
        },
        "review_status": review,
        "sequence": seq_block or {"value": seq, "length": len(seq) if seq else None},
        "isoforms": (up or {}).get("isoforms") or [],
        "sequence_conflicts": [],
        "sequence_analysis": seq_analysis,  # calculated
        "functional_annotations": (up or {}).get("functional_annotations") or {},
        "domains_motifs": (up or {}).get("domains_motifs") or [],
        "residue_annotations": (up or {}).get("residue_annotations") or [],
        "ptms": (up or {}).get("ptms") or [],
        "variants": (up or {}).get("variants") or [],
        "structures": {"experimental": [], "homologous": [], "predicted": [],
                       "best_for_docking": None, "assessed": False},
        "chemical_evidence": [],
        "protein_interactions": [],
        "pathways": [],
        "disease_associations": (up or {}).get("disease_associations") or [],
        "literature": [],
        "xref_graph": resolution.get("xref_graph") or {"nodes": [], "edges": []},
        "warnings": warnings,
        "provenance_summary": {
            "sources_used": sources_used,
            "retrieved_utc": _now(),
            "tool_version": f"SnaCleX v{__version__}",
            "identity_key": resolution.get("identity_key"),
        },
    }
    return record


def structure_availability(accession: str, *, sequence=None,
                           include_homologs: bool = False) -> dict:
    """Stage 4 hierarchy: experimental → homologous → predicted → sequence-only.

    Network-touching. Reports coverage/confidence and a docking-suitability gate;
    never auto-marks a low-confidence predicted model (or a distant homolog) as a
    reliable receptor. Homolog search is opt-in (``include_homologs`` + a
    ``sequence``) so the default availability check stays lightweight (NFR-4).
    """
    result = {"experimental": [], "homologous": [], "predicted": [],
              "tier": "sequence_only", "best_for_docking": None,
              "assessed": True, "warnings": []}

    # 1. Direct experimental structures via RCSB by-UniProt search.
    try:
        entities = rcsb.search_by_uniprot(accession)
    except FetchError as exc:
        entities = []
        result["warnings"].append(f"RCSB search unavailable: {exc}")
    seen = set()
    for e in entities:
        pid = e["pdb_id"]
        if pid in seen:
            continue
        seen.add(pid)
        result["experimental"].append({
            "kind": "experimental", "source": "RCSB", "pdb_id": pid,
            "entity": e.get("entity"), "docking_suitable": True,
            "provenance": {"source": "RCSB PDB", "source_id": pid,
                           "retrieved_utc": _now(), "evidence_category": "experimental"},
        })

    # 2. Homologous experimental structures (only when no direct structure and
    #    a sequence is available to search with).
    if include_homologs and sequence and not result["experimental"]:
        homo = homology.find_homologous_structures(
            sequence, direct_pdb_ids=[e["pdb_id"] for e in result["experimental"]])
        result["homologous"] = homo.get("structures", [])
        result["homology_method"] = homo.get("method")
        if not homo.get("available"):
            result["warnings"].append(
                f"homolog search unavailable: {homo.get('reason')}")

    # 3. Predicted structures via AlphaFold DB (metadata only — cheap).
    try:
        models = alphafold.fetch_models(accession)
    except FetchError as exc:
        models = []
        result["warnings"].append(f"AlphaFold unavailable: {exc}")
    for m in models:
        model = m["data"]
        model.update(alphafold.docking_assessment(model))
        model["provenance"] = m["provenance"]
        result["predicted"].append(model)

    # 4. Decide the tier (experimental > homologous > predicted > sequence-only).
    if result["experimental"]:
        result["tier"] = "experimental"
        result["best_for_docking"] = result["experimental"][0]["pdb_id"]
    elif result["homologous"]:
        result["tier"] = "homologous"
        dockable = next((h for h in result["homologous"] if h.get("docking_suitable")), None)
        if dockable:
            result["best_for_docking"] = dockable["pdb_id"]
            result["warnings"].append(
                "best available receptor is a HOMOLOG, not the requested protein — "
                "binding-site transfer requires caution")
        else:
            result["warnings"].append(
                "only distant homolog structures found — docking not offered")
    elif result["predicted"]:
        result["tier"] = "predicted"
        confident = next((m for m in result["predicted"] if m.get("docking_suitable")), None)
        if confident:
            result["best_for_docking"] = confident["model_id"]
        else:
            result["warnings"].append(
                "predicted model available but low confidence — docking not offered; "
                "sequence-only analysis remains available")
    else:
        result["warnings"].append(
            "no experimental, homologous, or predicted structure — sequence-only analysis")
    return result


def attach_structures(record: dict, availability: dict) -> dict:
    """Merge a structure-availability result into a record's structures block."""
    record["structures"] = {
        "experimental": availability.get("experimental", []),
        "homologous": availability.get("homologous", []),
        "predicted": availability.get("predicted", []),
        "best_for_docking": availability.get("best_for_docking"),
        "tier": availability.get("tier"),
        "assessed": True,
    }
    record.setdefault("warnings", []).extend(availability.get("warnings", []))
    if availability.get("tier") == "predicted":
        record["warnings"].append("structure is a predicted model — not experimental")
    return record
