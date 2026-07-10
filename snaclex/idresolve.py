"""Identifier resolution — Stage 1 (query interpretation) and Stage 2
(identifier normalization / cross-reference graph).

``interpret`` is pure and heavily unit-tested: it classifies any accepted input
without touching the network. ``resolve`` composes the adapters to build a
cross-reference graph and pick a candidate, degrading gracefully when a bridge is
unavailable. Identity is decided by accession + taxon + sequence checksum +
isoform role — never by gene symbol or name alone (see docs/platform/05).
"""

from __future__ import annotations

import re

from . import checksum, ncbi, rcsb, uniparc, uniprot
from .http_util import FetchError

# --- classification patterns ------------------------------------------------
_UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
    r"(?:-\d+)?$")
_UPI_RE = re.compile(r"^UPI[0-9A-F]{10}$")
_REFSEQ_RE = re.compile(r"^(?:NP|XP|WP|YP|AP)_\d+(?:\.\d+)?$", re.IGNORECASE)
_NCBI_RE = re.compile(r"^[A-Z]{3}\d{5,7}(?:\.\d+)?$")
_PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_VARIANT_RE = re.compile(r"^(?P<gene>[A-Za-z0-9]+)\s+(?P<wt>[A-Za-z])(?P<pos>\d+)(?P<mut>[A-Za-z*])$")
_HGVS_RE = re.compile(r"^(?P<acc>\S+):p\.[A-Za-z]{1,3}\d+[A-Za-z]{1,3}$")
_GENE_RE = re.compile(r"^[A-Z][A-Z0-9\-]{0,9}$")
_AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBZXUO*]+$", re.IGNORECASE)
_SMILES_HINT = set("=#()[]/\\@+.")


def _looks_chemical(token: str) -> bool:
    t = token.strip()
    if not t:
        return False
    if _INCHIKEY_RE.match(t) or t.startswith("InChI="):
        return True
    # SMILES: contains structural tokens, or organic-subset with lowercase atoms.
    if any(ch in _SMILES_HINT for ch in t) and not t[0].isspace():
        return True
    return False


def _classify_single(token: str) -> str:
    t = token.strip()
    if not t:
        return "unknown"
    if _UNIPROT_RE.match(t):
        return "uniprot_accession"
    if _UPI_RE.match(t):
        return "uniparc_upi"
    if _REFSEQ_RE.match(t):
        return "refseq_protein"
    if _INCHIKEY_RE.match(t) or t.startswith("InChI="):
        return "chemical"
    if _PDB_RE.match(t) and not _UNIPROT_RE.match(t):
        return "pdb_id"
    if _NCBI_RE.match(t):
        return "ncbi_protein"
    if _looks_chemical(t):
        return "chemical"
    if t.isdigit():
        return "chemical"  # PubChem CID in a chemical context
    if _GENE_RE.match(t):
        return "gene_symbol"
    return "protein_name"


def interpret(query: str) -> dict:
    """Classify an input into a query type with tokens (Stage 1). Pure/offline."""
    q = (query or "").strip()
    if not q:
        return {"query_type": "unknown", "tokens": [], "needs_disambiguation": False}

    # Raw FASTA (with header) or bare amino-acid sequence.
    if q.startswith(">"):
        return {"query_type": "raw_sequence", "tokens": [q], "needs_disambiguation": False}
    collapsed = "".join(q.split())
    if len(collapsed) >= 20 and _AA_RE.match(collapsed):
        return {"query_type": "raw_sequence", "tokens": [collapsed],
                "needs_disambiguation": False}

    # Variant: "BRAF V600E" or "P15056:p.Val600Glu".
    m = _VARIANT_RE.match(q)
    if m:
        return {"query_type": "variant", "tokens": [q],
                "variant": {"gene_or_acc": m.group("gene"), "wt": m.group("wt").upper(),
                            "position": int(m.group("pos")), "mut": m.group("mut").upper()},
                "needs_disambiguation": False}
    if _HGVS_RE.match(q):
        return {"query_type": "variant", "tokens": [q], "needs_disambiguation": False}

    # Explicit protein–chemical pair separated by '+'.
    if "+" in q:
        left, right = (s.strip() for s in q.split("+", 1))
        lt, rt = _classify_single(left), _classify_single(right)
        prot = {"uniprot_accession", "refseq_protein", "ncbi_protein", "pdb_id",
                "gene_symbol", "protein_name", "uniparc_upi"}
        if lt in prot and (rt == "chemical" or _looks_chemical(right)):
            return {"query_type": "protein_chemical_pair",
                    "tokens": [left, right],
                    "protein": {"query": left, "type": lt},
                    "chemical": {"query": right},
                    "needs_disambiguation": False}

    # Batch: multiple comma/newline separated identifiers.
    parts = [p.strip() for p in re.split(r"[\n,]", q) if p.strip()]
    if len(parts) > 1:
        return {"query_type": "batch", "tokens": parts,
                "item_types": [_classify_single(p) for p in parts],
                "needs_disambiguation": False}

    qtype = _classify_single(q)
    needs = qtype in ("gene_symbol", "protein_name")
    return {"query_type": qtype, "tokens": [q], "needs_disambiguation": needs}


# --- cross-reference graph --------------------------------------------------
def _node(ns, _id, qualifier=None):
    n = {"ns": ns, "id": _id}
    if qualifier:
        n["qualifier"] = qualifier
    return n


def build_xref_graph(*, uniprot_fragment=None, uniparc_fragment=None,
                     refseq_accessions=None, pdb_entities=None,
                     alphafold_ids=None) -> dict:
    """Assemble a cross-reference multigraph from parsed adapter fragments.

    Pure: takes already-fetched fragments so it is unit-testable offline. Keeps
    1:N and N:N links; never collapses isoforms/paralogs/orthologs.
    """
    nodes: list[dict] = []
    edges: list[dict] = []
    seen = set()

    def add_node(ns, _id, qualifier=None):
        key = (ns, _id, qualifier)
        if _id and key not in seen:
            seen.add(key)
            nodes.append(_node(ns, _id, qualifier))

    def add_edge(a, b, rel, source):
        edges.append({"from": a, "to": b, "rel": rel, "source": source})

    up_acc = None
    if uniprot_fragment:
        up_acc = uniprot_fragment.get("accession")
        add_node("UniProt", up_acc)
        xrefs = uniprot_fragment.get("cross_references") or {}
        for db, entries in xrefs.items():
            for e in entries:
                add_node(db, e["id"])
                add_edge(f"UniProt:{up_acc}", f"{db}:{e['id']}", "maps_to", "UniProtKB")

    if uniparc_fragment:
        upi = uniparc_fragment.get("upi")
        add_node("UniParc", upi)
        if up_acc:
            add_edge(f"UniProt:{up_acc}", f"UniParc:{upi}", "sequence_identical", "UniParc")
        for acc in uniparc.uniprot_accessions(uniparc_fragment):
            add_node("UniProt", acc)
            add_edge(f"UniParc:{upi}", f"UniProt:{acc}", "sequence_identical", "UniParc")

    for acc in refseq_accessions or []:
        add_node("RefSeq", acc)
        if up_acc:
            add_edge(f"UniProt:{up_acc}", f"RefSeq:{acc}", "maps_to", "resolver")

    for ent in pdb_entities or []:
        pdb, chain = (ent if isinstance(ent, (list, tuple)) else (ent, None))
        add_node("PDB", pdb, chain)
        if up_acc:
            add_edge(f"UniProt:{up_acc}", f"PDB:{pdb}", "has_structure", "resolver")

    for af in alphafold_ids or []:
        add_node("AlphaFold", af)
        if up_acc:
            add_edge(f"UniProt:{up_acc}", f"AlphaFold:{af}", "has_predicted_structure", "AlphaFold DB")

    return {"nodes": nodes, "edges": edges}


def _bridge_sequence_to_uniprot(seq: str) -> dict:
    """Raw sequence → CRC-64 → UniParc → UniProt (the checksum bridge)."""
    sums = checksum.checksums(seq)
    result = {"checksums": sums, "uniparc": None, "uniprot_accessions": []}
    try:
        hits = uniparc.search_by_checksum(sums["crc64"], limit=1)
    except FetchError:
        hits = []
    if hits:
        result["uniparc"] = hits[0]
        result["uniprot_accessions"] = uniparc.uniprot_accessions(hits[0])
    return result


def resolve(query: str, *, taxon=None, accession=None) -> dict:
    """Stage 2: normalize identifiers and build the cross-reference graph.

    Network-touching; degrades gracefully. Returns chosen entity + candidates +
    xref graph + identity key + warnings. Chemical inputs are returned with a
    flag so the caller routes them to the chemical path instead.
    """
    interp = interpret(query)
    qtype = interp["query_type"]
    warnings: list[str] = []
    out = {"query": query, "query_type": qtype, "chosen": None, "candidates": [],
           "xref_graph": {"nodes": [], "edges": []}, "identity_key": None,
           "warnings": warnings}

    if qtype == "chemical":
        out["is_chemical"] = True
        return out

    up_frag = None
    upi_frag = None
    refseq_accs: list[str] = []
    chosen_acc = accession

    try:
        if qtype in ("gene_symbol", "protein_name"):
            candidates = uniprot.search(query, organism_id=taxon, limit=10)
            out["candidates"] = candidates
            if not candidates:
                warnings.append("no UniProt candidates for this query")
                return out
            chosen_acc = accession or candidates[0]["accession"]
            if len(candidates) > 1:
                warnings.append(
                    f"{len(candidates)} candidates matched; showing "
                    f"{candidates[0].get('organism')} ({candidates[0].get('review_status')}) "
                    "— others available")
        elif qtype == "uniprot_accession":
            chosen_acc = query
        elif qtype == "uniparc_upi":
            upi_frag = uniparc.fetch_entry(query)["data"]
            accs = uniparc.uniprot_accessions(upi_frag)
            chosen_acc = accession or (accs[0] if accs else None)
        elif qtype in ("refseq_protein", "ncbi_protein"):
            nrec = ncbi.fetch_protein(query)["data"]
            refseq_accs = [nrec["accession"]] if nrec.get("accession") else []
            seq = (nrec.get("sequence") or {}).get("value")
            if seq:
                bridge = _bridge_sequence_to_uniprot(seq)
                out["identity_key"] = {"crc64": bridge["checksums"]["crc64"]}
                if bridge["uniparc"]:
                    upi_frag = bridge["uniparc"]
                accs = bridge["uniprot_accessions"]
                chosen_acc = accession or (accs[0] if accs else None)
                if not chosen_acc:
                    warnings.append("NCBI record has no UniProt mapping; sequence-only")
            out["ncbi"] = nrec
        elif qtype == "pdb_id":
            accs = rcsb.fetch_uniprot_accessions(query)
            chosen_acc = accession or (accs[0] if accs else None)
            if not chosen_acc:
                warnings.append("PDB entry has no UniProt reference")
        elif qtype == "raw_sequence":
            seq = _extract_sequence(query)
            bridge = _bridge_sequence_to_uniprot(seq)
            out["identity_key"] = {"crc64": bridge["checksums"]["crc64"]}
            out["sequence"] = seq
            if bridge["uniparc"]:
                upi_frag = bridge["uniparc"]
            accs = bridge["uniprot_accessions"]
            chosen_acc = accession or (accs[0] if accs else None)
            if not chosen_acc:
                warnings.append("novel/unmatched sequence — no UniProt/UniParc record; sequence-only")
    except FetchError as exc:
        warnings.append(f"resolution partially failed: {exc}")

    if chosen_acc:
        try:
            up_frag = uniprot.fetch_entry(chosen_acc)["data"]
        except FetchError as exc:
            warnings.append(f"could not fetch UniProt {chosen_acc}: {exc}")

    if up_frag:
        xrefs = up_frag.get("cross_references") or {}
        for e in xrefs.get("RefSeq", []):
            if e["id"] not in refseq_accs:
                refseq_accs.append(e["id"])
        pdb_entities = [(e["id"], None) for e in xrefs.get("PDB", [])]
        af_ids = [f"AF-{up_frag['accession']}-F1"] if up_frag.get("accession") else []
        if upi_frag is None:
            try:
                hits = uniparc.search_by_uniprot(up_frag["accession"], limit=1)
                upi_frag = hits[0] if hits else None
            except FetchError:
                pass
        out["xref_graph"] = build_xref_graph(
            uniprot_fragment=up_frag, uniparc_fragment=upi_frag,
            refseq_accessions=refseq_accs, pdb_entities=pdb_entities,
            alphafold_ids=af_ids)
        out["chosen"] = {
            "accession": up_frag.get("accession"),
            "taxon_id": (up_frag.get("organism") or {}).get("taxon_id"),
            "review_status": up_frag.get("review_status"),
        }
        out["identity_key"] = {
            "taxon_id": (up_frag.get("organism") or {}).get("taxon_id"),
            "crc64": (up_frag.get("sequence") or {}).get("crc64")
                     or (out.get("identity_key") or {}).get("crc64"),
            "isoform": "canonical",
        }
        out["_uniprot_fragment"] = up_frag
        out["_uniparc_fragment"] = upi_frag
    return out


def _extract_sequence(fasta_or_seq: str) -> str:
    """Return the amino-acid string from a FASTA record or a bare sequence."""
    lines = fasta_or_seq.splitlines()
    if lines and lines[0].startswith(">"):
        return "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    return "".join(fasta_or_seq.split())
