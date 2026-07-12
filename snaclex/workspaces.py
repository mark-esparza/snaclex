"""Research-domain workspace lenses (Phase 3).

Immunology / Oncology / Genetics *views* over the same ProteinRecord — they add
no new data source, they classify what the record already carries. Classification
is grounded first in **curated UniProt keywords** (authoritative), then in
documented gene/name/domain heuristics as a supplement. A small built-in seed
list of canonical examples is used only as a fallback and is explicitly labelled
as illustrative, not authoritative — absence from it never means the role is
absent.

Scientific safeguards enforced here:
* HLA/MHC allele identity is preserved and never collapsed under the gene symbol.
* The oncology lens separates "protein associated with cancer" from any claim
  about a compound being therapeutically effective — it makes neither a
  therapeutic nor a clinical recommendation.
"""

from __future__ import annotations

import re

# --- documented seed lists (illustrative canonical examples, NOT authoritative) ---
_ONCOGENE_SEED = {
    "MYC", "MYCN", "KRAS", "NRAS", "HRAS", "BRAF", "EGFR", "ERBB2", "ABL1",
    "ALK", "RET", "KIT", "MET", "PIK3CA", "MDM2", "CCND1", "CDK4", "JAK2",
    "FLT3", "SRC", "AKT1", "CTNNB1",
}
_TUMOR_SUPPRESSOR_SEED = {
    "TP53", "RB1", "PTEN", "APC", "BRCA1", "BRCA2", "VHL", "CDKN2A", "NF1",
    "NF2", "STK11", "SMAD4", "MLH1", "MSH2", "MSH6", "PMS2", "ATM", "WT1",
}
_DNA_REPAIR_SEED = {
    "BRCA1", "BRCA2", "ATM", "ATR", "MLH1", "MSH2", "MSH6", "PMS2", "RAD51",
    "PALB2", "XRCC1", "ERCC1", "FANCA",
}
_CHECKPOINT_SEED = {
    "PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "VSIR",
    "BTLA", "CD276", "CD80", "CD86",
}
_ANTIGEN_PROCESSING_SEED = {
    "TAP1", "TAP2", "TAPBP", "ERAP1", "ERAP2", "PSMB8", "PSMB9", "CALR", "CANX",
}

_CANCER_TERMS = re.compile(
    r"cancer|carcinoma|tumou?r|leuk[ae]mia|lymphoma|melanoma|sarcoma|glioma|"
    r"neoplasm|adenoma|blastoma|myeloma", re.IGNORECASE)
_HLA_ALLELE_RE = re.compile(r"\bHLA-[A-Z]+\d*\*\d+(?::\d+)*", re.IGNORECASE)


def _keywords(record: dict) -> set[str]:
    return {(k.get("name") or "").lower() for k in (record.get("keywords") or [])}


def _has_kw(kws: set[str], *targets) -> bool:
    return any(t.lower() in kws for t in targets)


def _gene(record: dict) -> str:
    return ((record.get("gene") or {}).get("symbol") or "").upper()


def _name(record: dict) -> str:
    return (record.get("preferred_name") or "").lower()


def _domains(record: dict) -> list[str]:
    return [(d.get("name") or "").lower() for d in (record.get("domains_motifs") or [])]


def _gene_prefix(gene: str, *prefixes) -> bool:
    return any(gene.startswith(p) for p in prefixes)


# ---------------------------------------------------------------------------
# Immunology
# ---------------------------------------------------------------------------
def immunology_view(record: dict, *, query: str | None = None) -> dict:
    kws = _keywords(record)
    gene = _gene(record)
    name = _name(record)
    categories: list[str] = []

    if _has_kw(kws, "cytokine") or _gene_prefix(gene, "IL", "IFN", "TNF", "CSF", "TGFB") \
            or "interleukin" in name or "interferon" in name:
        categories.append("cytokine / cytokine receptor")
    if _gene_prefix(gene, "CCL", "CXCL", "CX3CL", "XCL", "CCR", "CXCR", "CX3CR", "XCR") \
            or "chemokine" in name:
        categories.append("chemokine")
    if gene in _CHECKPOINT_SEED:
        categories.append("immune checkpoint")
    if _gene_prefix(gene, "IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD") \
            or "immunoglobulin" in name or any("immunoglobulin" in d for d in _domains(record)):
        categories.append("immunoglobulin / receptor")
    if gene in _ANTIGEN_PROCESSING_SEED:
        categories.append("antigen processing/presentation")
    if _has_kw(kws, "innate immunity", "adaptive immunity", "immunity",
               "inflammatory response"):
        categories.append("immune-related (keyword)")

    hla = _hla_mhc(record, query)
    if hla["is_hla_mhc"]:
        categories.append("HLA / MHC")

    return {
        "in_scope": bool(categories) or hla["is_hla_mhc"],
        "categories": sorted(set(categories)),
        "hla_mhc": hla,
        "immune_keywords": sorted(k for k in kws if any(
            t in k for t in ("immun", "cytokine", "inflammat", "mhc", "antigen"))),
        "small_molecules_note": ("Compounds affecting this protein are retrieved "
                                 "via /api/evidence (typed A–F), not asserted here."),
        "method": _method("immunology"),
    }


def _hla_mhc(record: dict, query: str | None) -> dict:
    gene = _gene(record)
    name = _name(record)
    kws = _keywords(record)
    is_hla = (gene.startswith("HLA-") or gene == "B2M"
              or "histocompatibility" in name or "mhc class" in name
              or _has_kw(kws, "mhc i", "mhc ii"))
    allele = None
    if query:
        m = _HLA_ALLELE_RE.search(query)
        if m:
            allele = m.group(0).upper()
    return {
        "is_hla_mhc": bool(is_hla),
        "allele_in_query": allele,
        "note": ("HLA/MHC gene: allele-level identifiers (e.g. HLA-A*02:01) are "
                 "functionally distinct and are NOT collapsed under the gene "
                 "symbol. This record is the gene-level canonical sequence; "
                 "allele-specific sequence handling requires the allele identifier."
                 ) if is_hla else None,
    }


# ---------------------------------------------------------------------------
# Oncology
# ---------------------------------------------------------------------------
def oncology_view(record: dict) -> dict:
    kws = _keywords(record)
    gene = _gene(record)
    name = _name(record)
    domains = _domains(record)
    roles: list[dict] = []

    def add(role, basis):
        roles.append({"role": role, "basis": basis})

    if _has_kw(kws, "proto-oncogene"):
        add("oncogene", "UniProt keyword 'Proto-oncogene' (curated)")
    elif gene in _ONCOGENE_SEED:
        add("oncogene", "canonical-example seed list (heuristic)")
    if _has_kw(kws, "tumor suppressor"):
        add("tumor suppressor", "UniProt keyword 'Tumor suppressor' (curated)")
    elif gene in _TUMOR_SUPPRESSOR_SEED:
        add("tumor suppressor", "canonical-example seed list (heuristic)")
    if _has_kw(kws, "kinase", "tyrosine-protein kinase", "serine/threonine-protein kinase") \
            or any("kinase" in d for d in domains) or "kinase" in name:
        add("kinase", "UniProt keyword/domain 'kinase'")
    if _has_kw(kws, "dna damage", "dna repair") or gene in _DNA_REPAIR_SEED:
        add("DNA repair", "UniProt keyword / seed list")
    if _has_kw(kws, "apoptosis"):
        add("apoptosis", "UniProt keyword 'Apoptosis' (curated)")
    if _has_kw(kws, "cell cycle"):
        add("cell cycle", "UniProt keyword 'Cell cycle' (curated)")
    if _has_kw(kws, "transcription", "transcription regulation", "activator", "repressor") \
            or any("dna-binding" in d or "zinc finger" in d for d in domains):
        add("transcription factor", "UniProt keyword/domain")

    cancer_diseases = [
        d for d in (record.get("disease_associations") or [])
        if _CANCER_TERMS.search((d.get("name") or "") + " " + (d.get("description") or ""))]
    driver_variants = [
        v for v in (record.get("variants") or [])
        if _CANCER_TERMS.search(v.get("description") or "")]

    return {
        "in_scope": bool(roles) or bool(cancer_diseases),
        "roles": roles,
        "cancer_associated": bool(cancer_diseases) or any(
            r["role"] in ("oncogene", "tumor suppressor") for r in roles),
        "cancer_disease_associations": [
            {"name": d.get("name"), "description": d.get("description"),
             "source": "UniProt (curated)"} for d in cancer_diseases],
        "cancer_associated_variants": [
            {"position": v.get("position"), "wt": v.get("wt"), "mut": v.get("mut"),
             "description": v.get("description")} for v in driver_variants],
        "structure_coverage_note": ("WT vs mutant structure coverage is available "
                                    "via /api/structure_availability and /api/variant."),
        "separation_note": ("This lens reports the protein's association with "
                            "cancer biology ONLY. It makes NO claim that any "
                            "compound is therapeutically effective; drug evidence "
                            "is typed separately in /api/evidence (levels A–F)."),
        "method": _method("oncology"),
    }


# ---------------------------------------------------------------------------
# Genetics (light summary; the residue-level path is /api/variant)
# ---------------------------------------------------------------------------
def genetics_view(record: dict) -> dict:
    variants = record.get("variants") or []
    return {
        "gene": (record.get("gene") or {}).get("symbol"),
        "gene_synonyms": (record.get("gene") or {}).get("synonyms") or [],
        "isoform_count": len(record.get("isoforms") or []),
        "variant_count": len(variants),
        "disease_associated_variants": [
            v for v in variants if v.get("description")][:50],
        "disease_associations": record.get("disease_associations") or [],
        "coding_consequence_note": ("Per-variant residue mapping + consequence via "
                                    "/api/variant. Frameshift/splice/start-loss/"
                                    "stop-loss and population frequency need "
                                    "transcript-level (Ensembl/VEP) and ClinVar."),
        "clinical_note": ("Variant interpretations are displayed with source and "
                          "review status only — never as a clinical recommendation."),
        "method": _method("genetics"),
    }


def _method(view: str) -> dict:
    return {
        "view": view,
        "basis": "curated UniProt keywords (primary) + documented gene/name/domain "
                 "heuristics + illustrative seed lists (supplementary)",
        "limitations": "Seed lists are canonical examples, not an authoritative "
                       "database; absence does not imply the role is absent. "
                       "Classifications are context tags, not clinical assertions.",
    }


def apply_views(record: dict, views, *, query: str | None = None) -> dict:
    out: dict = {}
    if "immunology" in views:
        out["immunology"] = immunology_view(record, query=query)
    if "oncology" in views:
        out["oncology"] = oncology_view(record)
    if "genetics" in views:
        out["genetics"] = genetics_view(record)
    return out
