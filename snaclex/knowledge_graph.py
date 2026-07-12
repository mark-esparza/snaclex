"""Provenance-aware knowledge graph assembly.

Projects a ProteinRecord (plus its structure availability, cross-reference graph,
and any typed evidence) into the node/edge schema from
docs/platform/04-data-model.md. The graph is a *projection* of already-resolved
data — it introduces no new claims — and **every edge carries provenance**
(source database, source id, evidence type, experimental-vs-predicted, retrieval
date, confidence, and software version for computed edges).

Pure and offline-testable: ``build_graph`` takes already-built inputs.
"""

from __future__ import annotations

from . import __version__

NODE_TYPES = [
    "Protein", "ProteinSequence", "Isoform", "Gene", "Organism", "Structure",
    "Domain", "Variant", "Compound", "Assay", "Disease", "Publication",
    "DatabaseRecord",
]
EDGE_TYPES = [
    "encoded_by", "in_organism", "has_sequence", "isoform_of",
    "sequence_identical_to", "homolog_of", "has_structure",
    "has_predicted_structure", "contains_domain", "has_variant", "binds",
    "inhibits", "activates", "interacts_with", "associated_with",
    "supported_by", "derived_from", "maps_to",
]

# Map an evidence predicate to an edge relation (all in EDGE_TYPES).
_PREDICATE_TO_REL = {
    "binds": "binds", "inhibits": "inhibits", "activates": "activates",
    "interacts_with": "interacts_with", "associated_with": "associated_with",
}


class _Graph:
    def __init__(self):
        self._nodes: dict[str, dict] = {}
        self.edges: list[dict] = []

    def node(self, node_id, ntype, label=None, **props):
        if node_id and node_id not in self._nodes:
            self._nodes[node_id] = {"id": node_id, "type": ntype,
                                    "label": label, "props": props}
        return node_id

    def edge(self, subj, rel, obj, *, source_database, source_id=None,
             evidence_type="curated", experimental_or_predicted="curated",
             publication=None, retrieval_date=None, confidence=None,
             software_version=None):
        if not subj or not obj:
            return
        self.edges.append({
            "from": subj, "to": obj, "rel": rel,
            "provenance": {
                "source_database": source_database,
                "source_id": source_id,
                "evidence_type": evidence_type,
                "experimental_or_predicted": experimental_or_predicted,
                "publication": publication,
                "retrieval_date": retrieval_date,
                "confidence_or_quality": confidence,
                "software_version": software_version,
            },
        })

    def result(self):
        nodes = list(self._nodes.values())
        node_counts: dict[str, int] = {}
        for n in nodes:
            node_counts[n["type"]] = node_counts.get(n["type"], 0) + 1
        edge_counts: dict[str, int] = {}
        for e in self.edges:
            edge_counts[e["rel"]] = edge_counts.get(e["rel"], 0) + 1
        return {
            "nodes": nodes, "edges": self.edges,
            "counts": {"nodes": len(nodes), "edges": len(self.edges),
                       "by_node_type": node_counts, "by_edge_type": edge_counts},
            "schema": {"node_types": NODE_TYPES, "edge_types": EDGE_TYPES},
            "note": ("Projection of resolved data; introduces no new claims. Every "
                     "edge carries provenance. Computed edges name their software."),
        }


def build_graph(record: dict, *, availability=None, evidence=None) -> dict:
    """Assemble the provenance-aware knowledge graph for a protein record."""
    g = _Graph()
    retrieved = ((record.get("provenance_summary") or {}).get("retrieved_utc"))
    accs = record.get("accessions") or {}
    protein_id = record.get("canonical_id")
    review = record.get("review_status")
    exp_status = "curated" if review in ("reviewed", "unreviewed") else review

    g.node(protein_id, "Protein", label=record.get("preferred_name"),
           accession=accs.get("uniprot_primary"), review_status=review,
           length=(record.get("sequence") or {}).get("length"))

    # Gene / organism.
    gene = (record.get("gene") or {}).get("symbol")
    if gene:
        gid = g.node(f"GENE:{gene}", "Gene", label=gene)
        g.edge(protein_id, "encoded_by", gid, source_database="UniProtKB",
               source_id=accs.get("uniprot_primary"), retrieval_date=retrieved)
    org = record.get("organism") or {}
    if org.get("taxon_id"):
        oid = g.node(f"TAXON:{org['taxon_id']}", "Organism",
                     label=org.get("scientific_name"))
        g.edge(protein_id, "in_organism", oid, source_database="UniProtKB",
               source_id=accs.get("uniprot_primary"), retrieval_date=retrieved)

    # Sequence node (identity).
    checksum = (record.get("sequence") or {}).get("checksum") or {}
    crc = checksum.get("crc64") or (record.get("sequence") or {}).get("crc64")
    if crc:
        sid = g.node(f"SEQ:{crc}", "ProteinSequence", crc64=crc,
                     length=(record.get("sequence") or {}).get("length"))
        g.edge(protein_id, "has_sequence", sid, source_database="UniParc",
               source_id=accs.get("uniparc"), retrieval_date=retrieved)

    # Isoforms.
    for iso in record.get("isoforms") or []:
        if iso.get("id"):
            nid = g.node(f"ISO:{iso['id']}", "Isoform", label=iso.get("name"))
            g.edge(nid, "isoform_of", protein_id, source_database="UniProtKB",
                   source_id=iso["id"], retrieval_date=retrieved)

    # Domains.
    for d in record.get("domains_motifs") or []:
        name = d.get("name")
        if not name:
            continue
        did = g.node(f"DOMAIN:{name}:{d.get('start')}-{d.get('end')}", "Domain",
                     label=name, start=d.get("start"), end=d.get("end"))
        g.edge(protein_id, "contains_domain", did,
               source_database=d.get("source") or "UniProtKB",
               source_id=accs.get("uniprot_primary"), retrieval_date=retrieved)

    # Variants.
    for v in record.get("variants") or []:
        pos = v.get("position")
        label = f"{v.get('wt') or ''}{pos}{v.get('mut') or ''}"
        vid = g.node(f"VAR:{protein_id}:{label}", "Variant", label=label,
                     description=v.get("description"))
        g.edge(protein_id, "has_variant", vid, source_database="UniProtKB",
               source_id=accs.get("uniprot_primary"), retrieval_date=retrieved)

    # Disease associations.
    for dis in record.get("disease_associations") or []:
        if dis.get("name"):
            nid = g.node(f"DISEASE:{dis['name']}", "Disease", label=dis["name"])
            g.edge(protein_id, "associated_with", nid, source_database="UniProtKB",
                   source_id=dis.get("accession"), retrieval_date=retrieved)

    # Literature.
    for lit in record.get("literature") or []:
        ref = lit.get("pmid") or lit.get("doi")
        if not ref:
            continue
        pid = g.node(f"PMID:{lit['pmid']}" if lit.get("pmid") else f"DOI:{lit['doi']}",
                     "Publication", label=lit.get("title"))
        g.edge(protein_id, "supported_by", pid,
               source_database=lit.get("source") or "UniProtKB",
               publication=(f"PMID:{lit['pmid']}" if lit.get("pmid") else None),
               retrieval_date=retrieved)

    _add_structures(g, protein_id, record, availability, accs, retrieved)
    _add_xrefs(g, protein_id, record)
    _add_evidence(g, protein_id, evidence, retrieved)
    return g.result()


def _add_structures(g, protein_id, record, availability, accs, retrieved):
    structures = availability or record.get("structures") or {}
    for s in structures.get("experimental") or []:
        pid = s.get("pdb_id")
        if pid:
            nid = g.node(f"PDB:{pid}", "Structure", label=pid, kind="experimental")
            g.edge(protein_id, "has_structure", nid, source_database="RCSB PDB",
                   source_id=pid, evidence_type="experimental",
                   experimental_or_predicted="experimental", retrieval_date=retrieved)
    for s in structures.get("homologous") or []:
        pid = s.get("pdb_id")
        if pid:
            nid = g.node(f"PDB:{pid}", "Structure", label=pid, kind="homologous",
                         sequence_identity=s.get("sequence_identity"))
            g.edge(protein_id, "has_structure", nid, source_database="RCSB PDB",
                   source_id=pid, evidence_type="homology-transferred",
                   experimental_or_predicted="experimental (homolog)",
                   confidence=s.get("sequence_identity"), retrieval_date=retrieved)
    for m in structures.get("predicted") or []:
        mid = m.get("model_id")
        if mid:
            nid = g.node(f"AF:{mid}", "Structure", label=mid, kind="predicted",
                         mean_plddt=m.get("global_plddt_mean"))
            g.edge(protein_id, "has_predicted_structure", nid,
                   source_database="AlphaFold DB", source_id=mid,
                   evidence_type="predicted", experimental_or_predicted="predicted",
                   confidence=m.get("global_plddt_mean"), retrieval_date=retrieved)


def _add_xrefs(g, protein_id, record):
    xg = record.get("xref_graph") or {}
    for n in xg.get("nodes") or []:
        ns, xid = n.get("ns"), n.get("id")
        if ns and xid and ns != "UniProt":
            g.node(f"{ns}:{xid}", "DatabaseRecord", label=f"{ns}:{xid}", namespace=ns)
    for e in xg.get("edges") or []:
        frm, to, rel = e.get("from"), e.get("to"), e.get("rel")
        # xref graph ids look like "UniProt:P00519"; map UniProt→the protein node.
        frm_id = protein_id if str(frm).startswith("UniProt:") else frm
        to_id = protein_id if str(to).startswith("UniProt:") else to
        mapped_rel = rel if rel in EDGE_TYPES else "maps_to"
        g.edge(frm_id, mapped_rel, to_id, source_database=e.get("source") or "resolver",
               evidence_type="cross-reference")


def _add_evidence(g, protein_id, evidence, retrieved):
    for ev in evidence or []:
        obj = ev.get("object")
        obj_ref = obj.get("ref") if isinstance(obj, dict) else obj
        if not obj_ref:
            continue
        cid = g.node(obj_ref, "Compound", label=(ev.get("chemical_form") or {}).get("inchikey"))
        rel = _PREDICATE_TO_REL.get(ev.get("predicate"), "binds")
        comp = ev.get("computation") or {}
        g.edge(protein_id, rel, cid,
               source_database=ev.get("source"),
               source_id=ev.get("source_record"),
               evidence_type=f"level {ev.get('level')} ({ev.get('category')})",
               experimental_or_predicted=("experimental" if ev.get("is_direct")
                                          else "predicted/inferred"),
               publication=ev.get("publication"),
               confidence=ev.get("confidence"),
               software_version=comp.get("software"),
               retrieval_date=ev.get("retrieved_utc") or retrieved)
        if ev.get("publication"):
            pid = g.node(ev["publication"], "Publication")
            g.edge(protein_id, "supported_by", pid,
                   source_database=ev.get("source"), publication=ev.get("publication"))
