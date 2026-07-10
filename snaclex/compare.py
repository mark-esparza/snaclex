"""Protein comparison — family / ortholog views over the shared record model.

Pure functions that take a list of already-built ProteinRecords and produce a
side-by-side comparison plus domain-sharing and orthology grouping. Used by the
batch job (``run_batch_job``) and any comparison surface.

Identity safeguard: proteins are **grouped**, never merged, and orthology is
proposed only as *candidates* keyed on shared gene symbol across differing taxa —
explicitly flagged as needing sequence/phylogeny confirmation, never treated as
established on symbol alone (per the task's identity rules).
"""

from __future__ import annotations

_COLUMNS = [
    "accession", "preferred_name", "gene", "organism", "taxon_id",
    "review_status", "length", "molecular_weight_Da", "isoelectric_point",
    "gravy", "n_domains", "structure_tier", "n_variants",
]


def record_summary(record: dict) -> dict:
    """One comparable row from a ProteinRecord."""
    sa = record.get("sequence_analysis") or {}
    gene = record.get("gene") or {}
    org = record.get("organism") or {}
    structures = record.get("structures") or {}
    domains = record.get("domains_motifs") or []
    return {
        "canonical_id": record.get("canonical_id"),
        "accession": (record.get("accessions") or {}).get("uniprot_primary"),
        "preferred_name": record.get("preferred_name"),
        "gene": gene.get("symbol"),
        "organism": org.get("scientific_name"),
        "taxon_id": org.get("taxon_id"),
        "review_status": record.get("review_status"),
        "length": (record.get("sequence") or {}).get("length") or sa.get("length"),
        "molecular_weight_Da": sa.get("molecular_weight_Da"),
        "isoelectric_point": sa.get("isoelectric_point"),
        "gravy": sa.get("gravy"),
        "n_domains": len(domains),
        "domain_names": sorted({d.get("name") for d in domains if d.get("name")}),
        "structure_tier": structures.get("tier"),
        "n_variants": len(record.get("variants") or []),
    }


def shared_domains(records: list[dict]) -> dict:
    """Domain names present in all records vs only some (family signal)."""
    sets = []
    for r in records:
        names = {d.get("name") for d in (r.get("domains_motifs") or []) if d.get("name")}
        sets.append(names)
    if not sets:
        return {"in_all": [], "in_some": []}
    in_all = set.intersection(*sets) if sets else set()
    in_any = set.union(*sets) if sets else set()
    return {"in_all": sorted(in_all), "in_some": sorted(in_any - in_all)}


def orthology_groups(records: list[dict]) -> dict:
    """Group by shared gene symbol → ortholog / paralog *candidates* (not merges)."""
    groups: dict[str, list] = {}
    for r in records:
        gene = (r.get("gene") or {}).get("symbol")
        if not gene:
            continue
        groups.setdefault(gene.upper(), []).append({
            "accession": (r.get("accessions") or {}).get("uniprot_primary"),
            "organism": (r.get("organism") or {}).get("scientific_name"),
            "taxon_id": (r.get("organism") or {}).get("taxon_id"),
        })
    orthologs = {g: m for g, m in groups.items()
                 if len({x["taxon_id"] for x in m}) > 1}
    paralogs = {g: m for g, m in groups.items()
                if len(m) > 1 and len({x["taxon_id"] for x in m}) == 1}
    return {
        "ortholog_candidates": orthologs,
        "paralog_candidates": paralogs,
        "note": ("Grouped by shared gene symbol only: ortholog = same symbol "
                 "across different organisms; paralog candidate = same symbol, "
                 "same organism, multiple records. CANDIDATES — confirm by "
                 "sequence identity / phylogeny, never by symbol alone."),
    }


def build_comparison(records: list[dict]) -> dict:
    """Full comparison payload over a list of ProteinRecords (pure)."""
    rows = [record_summary(r) for r in records]
    lengths = [row["length"] for row in rows if row["length"]]
    return {
        "n": len(rows),
        "columns": _COLUMNS,
        "rows": rows,
        "shared_domains": shared_domains(records),
        "orthology": orthology_groups(records),
        "length_range": ({"min": min(lengths), "max": max(lengths)}
                         if lengths else None),
        "note": ("Comparison over the shared ProteinRecord model. Sequence metrics "
                 "are calculated; identity/orthology grouping is by gene symbol and "
                 "is a candidate view, not a merge."),
    }
