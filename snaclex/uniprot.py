"""UniProtKB adapter — the primary protein-annotation layer.

Network entry points are thin wrappers around the documented REST API
(``https://rest.uniprot.org/uniprotkb/...``); the parsing is isolated in
``parse_entry``/``parse_search`` so it is unit-testable offline against recorded
fixtures (matching the repo's existing adapter-test convention).

Adapters return **normalized fragments** — source-agnostic keys plus a provenance
envelope — never raw UniProt JSON. Swiss-Prot (reviewed) vs TrEMBL (unreviewed)
is always surfaced. Fields absent upstream become ``None`` with a note; nothing
is fabricated. Paths marked ``VERIFY`` need a live round-trip before production
reliance (see docs/platform/03-source-integration.md).
"""

from __future__ import annotations

import datetime
import urllib.parse

from .http_util import FetchError, fetch_json

_BASE = "https://rest.uniprot.org/uniprotkb"

# Reviewed vs unreviewed as UniProt reports it in `entryType`.
_REVIEWED = "UniProtKB reviewed (Swiss-Prot)"
_UNREVIEWED = "UniProtKB unreviewed (TrEMBL)"


def _now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _provenance(acc: str, version=None, category="curated", limitations=None) -> dict:
    return {
        "source": "UniProtKB",
        "source_id": acc,
        "source_version": version,
        "retrieved_utc": _now(),
        "evidence_category": category,
        "limitations": limitations,
    }


# ---------------------------------------------------------------------------
# Network entry points (thin)
# ---------------------------------------------------------------------------
def fetch_entry(accession: str) -> dict:
    """Fetch and normalize a single UniProtKB entry by accession."""
    acc = urllib.parse.quote((accession or "").strip())
    if not acc:
        raise FetchError("empty UniProt accession")
    data = fetch_json(f"{_BASE}/{acc}.json")
    return parse_entry(data)


def search(query: str, *, organism_id=None, limit: int = 10) -> list[dict]:
    """Search UniProtKB, returning ranked candidate proteins (reviewed first)."""
    q = query
    if organism_id:
        q = f"({query}) AND organism_id:{organism_id}"
    fields = "accession,id,protein_name,gene_names,organism_name,organism_id,length,reviewed"
    url = (
        f"{_BASE}/search?query={urllib.parse.quote(q)}"
        f"&fields={fields}&format=json&size={int(limit)}"
    )
    data = fetch_json(url)
    return parse_search(data)


# ---------------------------------------------------------------------------
# Pure parsing (offline-testable)
# ---------------------------------------------------------------------------
def review_status(entry_type: str | None) -> str:
    if entry_type == _REVIEWED:
        return "reviewed"
    if entry_type == _UNREVIEWED:
        return "unreviewed"
    return "unknown"


def parse_search(data: dict) -> list[dict]:
    out: list[dict] = []
    for r in (data or {}).get("results", []):
        desc = r.get("proteinDescription") or {}
        rec = ((desc.get("recommendedName") or {}).get("fullName") or {}).get("value")
        genes = r.get("genes") or []
        gene = None
        if genes and genes[0].get("geneName"):
            gene = genes[0]["geneName"].get("value")
        org = r.get("organism") or {}
        out.append({
            "accession": r.get("primaryAccession"),
            "uniprot_id": r.get("uniProtkbId"),
            "protein_name": rec,
            "gene": gene,
            "organism": org.get("scientificName"),
            "taxon_id": org.get("taxonId"),
            "length": (r.get("sequence") or {}).get("length"),
            "review_status": review_status(r.get("entryType")),
        })
    return out


def _gene(data: dict) -> dict:
    genes = data.get("genes") or []
    if not genes:
        return {"symbol": None, "synonyms": []}
    g = genes[0]
    symbol = (g.get("geneName") or {}).get("value")
    synonyms = [s.get("value") for s in (g.get("synonyms") or []) if s.get("value")]
    return {"symbol": symbol, "synonyms": synonyms}


def _preferred_name(data: dict) -> str | None:
    desc = data.get("proteinDescription") or {}
    rec = (desc.get("recommendedName") or {}).get("fullName") or {}
    if rec.get("value"):
        return rec["value"]
    sub = desc.get("submissionNames") or []
    if sub:
        return ((sub[0].get("fullName") or {}).get("value"))
    return None


def _loc(feature: dict) -> tuple:
    loc = feature.get("location") or {}
    start = (loc.get("start") or {}).get("value")
    end = (loc.get("end") or {}).get("value")
    return start, end


_RESIDUE_KINDS = {
    "Signal": "signal_peptide",
    "Transit peptide": "transit_peptide",
    "Transmembrane": "transmembrane",
    "Intramembrane": "intramembrane",
    "Coiled coil": "coiled_coil",
    "Active site": "active_site",
    "Binding site": "binding_site",
    "Site": "site",
    "Metal binding": "binding_site",
    "Disulfide bond": "disulfide",
}
_PTM_KINDS = {"Modified residue", "Glycosylation", "Lipidation", "Cross-link"}
_DOMAIN_KINDS = {"Domain", "Region", "Repeat", "Motif", "Zinc finger", "DNA binding"}


def classify_features(data: dict) -> dict:
    """Split UniProt ``features`` into domains, residue annotations, PTMs, variants."""
    domains, residues, ptms, variants = [], [], [], []
    for f in data.get("features") or []:
        ftype = f.get("type")
        start, end = _loc(f)
        desc = f.get("description")
        if ftype in _DOMAIN_KINDS:
            domains.append({"type": ftype.lower(), "name": desc, "start": start, "end": end})
        elif ftype in _RESIDUE_KINDS:
            residues.append({"kind": _RESIDUE_KINDS[ftype], "start": start, "end": end,
                             "description": desc})
        elif ftype in _PTM_KINDS:
            ptms.append({"type": ftype.lower(), "position": start, "description": desc})
        elif ftype == "Natural variant":
            alt = f.get("alternativeSequence") or {}
            orig = alt.get("originalSequence")
            alts = alt.get("alternativeSequences") or []
            variants.append({
                "position": start, "wt": orig,
                "mut": alts[0] if alts else None,
                "kind": "missense" if orig and alts else "variant",
                "description": desc,
            })
    return {"domains": domains, "residues": residues, "ptms": ptms, "variants": variants}


def _comments(data: dict) -> dict:
    functions, catalytic, cofactors, localization, tissue = [], [], [], [], []
    isoforms, diseases = [], []
    for c in data.get("comments") or []:
        ctype = c.get("commentType")
        if ctype == "FUNCTION":
            functions += [t.get("value") for t in c.get("texts") or [] if t.get("value")]
        elif ctype == "CATALYTIC ACTIVITY":
            name = (c.get("reaction") or {}).get("name")
            if name:
                catalytic.append(name)
        elif ctype == "COFACTOR":
            cofactors += [cf.get("name") for cf in c.get("cofactors") or [] if cf.get("name")]
        elif ctype == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations") or []:
                val = (sl.get("location") or {}).get("value")
                if val:
                    localization.append(val)
        elif ctype == "TISSUE SPECIFICITY":
            tissue += [t.get("value") for t in c.get("texts") or [] if t.get("value")]
        elif ctype == "ALTERNATIVE PRODUCTS":
            for iso in c.get("isoforms") or []:
                ids = iso.get("isoformIds") or []
                isoforms.append({
                    "id": ids[0] if ids else None,
                    "name": (iso.get("name") or {}).get("value"),
                    "sequence_status": iso.get("isoformSequenceStatus"),
                })
        elif ctype == "DISEASE":
            d = c.get("disease") or {}
            if d.get("diseaseId") or d.get("diseaseAccession"):
                diseases.append({
                    "name": d.get("diseaseId"),
                    "acronym": d.get("acronym"),
                    "accession": d.get("diseaseAccession"),
                    "description": d.get("description"),
                })
    return {"functions": functions, "catalytic_activity": catalytic,
            "cofactors": cofactors, "subcellular_location": localization,
            "tissue_specificity": tissue, "isoforms": isoforms, "diseases": diseases}


def cross_references(data: dict) -> dict:
    """Collect useful cross-references keyed by database (RefSeq, PDB, Ensembl…)."""
    xrefs: dict[str, list] = {}
    for x in data.get("uniProtKBCrossReferences") or []:
        db = x.get("database")
        xid = x.get("id")
        if not db or not xid:
            continue
        xrefs.setdefault(db, [])
        if xid not in [e["id"] for e in xrefs[db]]:
            props = {p.get("key"): p.get("value") for p in x.get("properties") or []}
            xrefs[db].append({"id": xid, "properties": props})
    return xrefs


def parse_entry(data: dict) -> dict:
    """Normalize a UniProtKB entry JSON into a fragment + provenance."""
    acc = data.get("primaryAccession")
    seq = data.get("sequence") or {}
    org = data.get("organism") or {}
    audit = data.get("entryAudit") or {}
    status = review_status(data.get("entryType"))
    features = classify_features(data)
    comments = _comments(data)

    fragment = {
        "accession": acc,
        "secondary_accessions": data.get("secondaryAccessions") or [],
        "uniprot_id": data.get("uniProtkbId"),
        "review_status": status,
        "preferred_name": _preferred_name(data),
        "gene": _gene(data),
        "organism": {"scientific_name": org.get("scientificName"),
                     "taxon_id": org.get("taxonId")},
        "sequence": {
            "value": seq.get("value"),
            "length": seq.get("length"),
            "molecular_weight_Da": seq.get("molWeight"),
            "crc64": seq.get("crc64"),
            "md5": seq.get("md5"),
        },
        "functional_annotations": {
            "functions": comments["functions"],
            "catalytic_activity": comments["catalytic_activity"],
            "cofactors": comments["cofactors"],
            "subcellular_location": comments["subcellular_location"],
            "tissue_specificity": comments["tissue_specificity"],
        },
        "domains_motifs": features["domains"],
        "residue_annotations": features["residues"],
        "ptms": features["ptms"],
        "variants": features["variants"],
        "isoforms": comments["isoforms"],
        "disease_associations": comments["diseases"],
        "cross_references": cross_references(data),
    }
    return {
        "data": fragment,
        "provenance": _provenance(
            acc, version=audit.get("entryVersion"),
            limitations=None if status == "reviewed" else "TrEMBL, unreviewed"),
        "status": "ok",
    }
