"""RCSB Protein Data Bank client: structure files + entry metadata."""

from __future__ import annotations

import re

from .http_util import FetchError, fetch_json, fetch_text

_PDB_ID_RE = re.compile(r"^[0-9A-Za-z]{4}$")


def normalize_pdb_id(pdb_id: str) -> str:
    pid = (pdb_id or "").strip().upper()
    if not _PDB_ID_RE.match(pid):
        raise FetchError(f"'{pdb_id}' is not a valid 4-character PDB ID")
    return pid


def fetch_structure(pdb_id: str) -> str:
    """Return the raw PDB-format text for an entry."""
    pid = normalize_pdb_id(pdb_id)
    return fetch_text(f"https://files.rcsb.org/download/{pid}.pdb")


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
