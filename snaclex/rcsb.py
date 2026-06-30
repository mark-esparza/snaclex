"""RCSB Protein Data Bank client: structure files + entry metadata."""

from __future__ import annotations

import re

from .http_util import FetchError, RateLimitError, fetch_json, fetch_text

_PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")


def normalize_pdb_id(pdb_id: str) -> str:
    pid = (pdb_id or "").strip().upper()
    if not _PDB_ID_RE.match(pid):
        raise FetchError(f"'{pdb_id}' is not a valid 4-character PDB ID")
    return pid


def fetch_structure(pdb_id: str) -> str:
    """Return the raw structure text for an entry (PDB format, or mmCIF fallback).

    Many large/newer entries have no legacy PDB-format file and 404 on the
    ``.pdb`` download; for those RCSB serves only mmCIF, so we fall back to it.
    ``pdbparse.parse_structure`` auto-detects which format it got.
    """
    pid = normalize_pdb_id(pdb_id)
    try:
        return fetch_text(f"https://files.rcsb.org/download/{pid}.pdb")
    except RateLimitError:
        raise
    except FetchError:
        return fetch_text(f"https://files.rcsb.org/download/{pid}.cif")


def fetch_entry_metadata(pdb_id: str) -> dict:
    """Return a compact, UI-friendly metadata dict for a PDB entry."""
    pid = normalize_pdb_id(pdb_id)
    data = fetch_json(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")

    struct = data.get("struct") or {}
    entry_info = data.get("rcsb_entry_info") or {}
    accession = data.get("rcsb_accession_info") or {}
    exptl = data.get("exptl") or [{}]

    resolution = None
    res_list = entry_info.get("resolution_combined")
    if isinstance(res_list, list) and res_list:
        resolution = res_list[0]

    methods = [e.get("method") for e in exptl if e.get("method")]

    return {
        "pdb_id": pid,
        "title": struct.get("title"),
        "experimental_method": ", ".join(methods) if methods else None,
        "resolution_A": resolution,
        "deposited": accession.get("initial_release_date"),
        "polymer_entity_count": entry_info.get("polymer_entity_count"),
        "deposited_atom_count": entry_info.get("deposited_atom_count"),
        "deposited_model_count": entry_info.get("deposited_model_count"),
        "molecular_weight_kDa": entry_info.get("molecular_weight"),
        "nonpolymer_count": entry_info.get("nonpolymer_entity_count"),
    }


def search_by_name(query: str, limit: int = 10) -> list[dict]:
    """Full-text search the PDB, returning [{id, score}] ranked entries."""
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": limit}},
    }
    import json
    import urllib.parse

    url = (
        "https://search.rcsb.org/rcsbsearch/v2/query?json="
        + urllib.parse.quote(json.dumps(payload))
    )
    data = fetch_json(url)
    results = []
    for item in data.get("result_set", []):
        results.append({"pdb_id": item.get("identifier"), "score": item.get("score")})

    # Enrich with title + organism in a single batched GraphQL call.
    try:
        summaries = fetch_entry_summaries([r["pdb_id"] for r in results])
        for r in results:
            s = summaries.get(r["pdb_id"], {})
            r["title"] = s.get("title")
            r["organism"] = s.get("organism")
    except FetchError:
        pass  # search still works without the enrichment
    return results


def fetch_uniprot_accessions(pdb_id: str) -> list[str]:
    """Return the UniProt accessions referenced by a PDB entry (may be empty)."""
    pid = normalize_pdb_id(pdb_id)
    import urllib.parse

    query = (
        '{entry(entry_id:"' + pid + '"){polymer_entities{'
        "rcsb_polymer_entity_container_identifiers{"
        "reference_sequence_identifiers{database_name database_accession}}}}}"
    )
    url = "https://data.rcsb.org/graphql?query=" + urllib.parse.quote(query)
    try:
        data = fetch_json(url)
    except FetchError:
        return []
    accs: list[str] = []
    entry = (data.get("data") or {}).get("entry") or {}
    for pe in entry.get("polymer_entities") or []:
        ids = pe.get("rcsb_polymer_entity_container_identifiers") or {}
        for ref in ids.get("reference_sequence_identifiers") or []:
            if ref.get("database_name") == "UniProt" and ref.get("database_accession"):
                acc = ref["database_accession"]
                if acc not in accs:
                    accs.append(acc)
    return accs


def _parse_chemcomp_descriptors(descriptors) -> str | None:
    """First SMILES descriptor from a pdbx_chem_comp_descriptor list."""
    for d in descriptors or []:
        if "SMILES" in (d.get("type") or "").upper() and d.get("descriptor"):
            return d["descriptor"]
    return None


def _fetch_chem_component_rest(code: str) -> dict | None:
    """Resolve a single component via the stable REST chemcomp endpoint.

    Used as a fallback when the batched GraphQL query doesn't resolve a code,
    so names keep working even if the GraphQL schema drifts.
    """
    try:
        data = fetch_json(f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}")
    except FetchError:
        return None
    meta = data.get("chem_comp") or {}
    if not meta.get("name") and not meta.get("formula"):
        return None
    synonyms = [
        s.get("name") for s in (data.get("rcsb_chem_comp_synonyms") or [])
        if s.get("name")
    ]
    return {
        "name": meta.get("name"),
        "formula": meta.get("formula"),
        "formula_weight": meta.get("formula_weight"),
        "type": meta.get("type"),
        "smiles": _parse_chemcomp_descriptors(data.get("pdbx_chem_comp_descriptor")),
        "synonyms": synonyms[:5],
    }


def fetch_chem_components(comp_ids: list[str]) -> dict:
    """Resolve PDB chemical-component codes to names + chemistry via the CCD.

    Returns ``{CODE: {name, formula, formula_weight, smiles, synonyms, type}}``
    for the codes that resolve. Tries one batched RCSB GraphQL call (mirrors
    ``fetch_entry_summaries``), then falls back to the stable per-code REST
    endpoint for any codes GraphQL didn't resolve. Returns only what resolved,
    so a schema hiccup or offline run never blocks structure loading.
    """
    codes = sorted({(c or "").strip().upper() for c in comp_ids if (c or "").strip()})
    if not codes:
        return {}
    import urllib.parse

    out: dict[str, dict] = {}
    id_list = ",".join(f'"{c}"' for c in codes)
    query = (
        "{chem_comps(comp_ids:[" + id_list + "])"
        "{chem_comp{id name formula formula_weight type pdbx_synonyms}"
        "pdbx_chem_comp_descriptor{type descriptor}}}"
    )
    url = "https://data.rcsb.org/graphql?query=" + urllib.parse.quote(query)
    try:
        data = fetch_json(url)
    except FetchError:
        data = {}

    for cc in (data.get("data", {}).get("chem_comps") or []):
        meta = cc.get("chem_comp") or {}
        code = (meta.get("id") or "").upper()
        if not code:
            continue
        synonyms = [
            s.strip() for s in (meta.get("pdbx_synonyms") or "").split(";") if s.strip()
        ]
        out[code] = {
            "name": meta.get("name"),
            "formula": meta.get("formula"),
            "formula_weight": meta.get("formula_weight"),
            "type": meta.get("type"),
            "smiles": _parse_chemcomp_descriptors(cc.get("pdbx_chem_comp_descriptor")),
            "synonyms": synonyms[:5],
        }

    # REST fallback for anything GraphQL missed (schema drift, partial result).
    for code in codes:
        if code not in out or not out[code].get("name"):
            rest = _fetch_chem_component_rest(code)
            if rest:
                out[code] = rest
    return out


def fetch_entry_summaries(ids: list[str]) -> dict:
    """Batch-fetch {pdb_id: {title, organism}} for several entries at once."""
    if not ids:
        return {}
    import urllib.parse

    id_list = ",".join(f'"{i}"' for i in ids)
    query = (
        "{entries(entry_ids:[" + id_list + "])"
        "{rcsb_id struct{title} "
        "polymer_entities{rcsb_entity_source_organism{ncbi_scientific_name}}}}"
    )
    url = "https://data.rcsb.org/graphql?query=" + urllib.parse.quote(query)
    data = fetch_json(url)

    out: dict[str, dict] = {}
    for e in (data.get("data", {}).get("entries") or []):
        rid = e.get("rcsb_id")
        title = (e.get("struct") or {}).get("title")
        organism = None
        for pe in e.get("polymer_entities") or []:
            srcs = pe.get("rcsb_entity_source_organism") or []
            if srcs and srcs[0].get("ncbi_scientific_name"):
                organism = srcs[0]["ncbi_scientific_name"]
                break
        out[rid] = {"title": title, "organism": organism}
    return out
