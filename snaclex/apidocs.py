"""Machine-readable API contract for SnaCleX, served at GET /api/docs.

Keeping the contract in one place (rather than scattered in handler docstrings)
lets the UI render a reference page and lets tests assert every routed endpoint
is documented. This is a hand-written contract, not generated — SnaCleX is
dependency-free, so there's no framework to introspect.
"""

from __future__ import annotations

from . import __version__

# Limits mirror the constants enforced in server.py.
_LIMITS = {
    "max_query_length": 200,
    "max_batch_chems_length": 2000,
    "max_batch_chems": 10,
    "max_upload_bytes": 5_000_000,
    "rate_limit": "per-IP token bucket; 429 + Retry-After when exceeded",
    "job_concurrency": "bounded worker pool; heavy work queues via /api/jobs",
}

_ENDPOINTS = [
    {
        "method": "GET", "path": "/api/analyze",
        "params": {"pdb": "4-char PDB ID or upload id (required)",
                   "chain": "optional chain id to load just that chain"},
        "returns": "id, metadata, chains, components, protein_atom_count, "
                   "pdb_data. For assemblies over the atom limit (no chain "
                   "given): {too_large, n_atoms, limit, chains:[{chain, "
                   "atom_count}]} so a chain can be chosen.",
    },
    {
        "method": "GET", "path": "/api/interactions",
        "params": {"pdb": "structure id (required)", "comp": "component index (required)"},
        "returns": "atomic interaction profile + plain-language report",
    },
    {
        "method": "GET", "path": "/api/chemical",
        "params": {"q": "drug/chemical/element name or CID (required)",
                   "pdb": "optional structure id for ChEMBL target cross-reference"},
        "returns": "PubChem properties, Lipinski + druglikeness, ChEMBL status",
    },
    {
        "method": "GET", "path": "/api/pockets",
        "params": {"pdb": "structure id (required)"},
        "returns": "ranked geometric cavities + methods/provenance",
    },
    {
        "method": "GET", "path": "/api/evolution",
        "params": {"pdb": "structure id (required)"},
        "returns": "Pfam conservation per residue/pocket + methods/provenance",
    },
    {
        "method": "GET", "path": "/api/search",
        "params": {"q": "free-text query (required)"},
        "returns": "ranked PDB full-text search results",
    },
    {
        "method": "GET", "path": "/api/resolve",
        "params": {"q": "protein/gene/accession/PDB/FASTA/variant (required)",
                   "taxon": "optional NCBI taxon id filter",
                   "accession": "optional explicit accession to disambiguate"},
        "returns": "Stage 1/2: query_type, candidates, cross-reference graph, "
                   "identity key (accession+taxon+checksum). Sequence-first.",
    },
    {
        "method": "GET", "path": "/api/protein",
        "params": {"q": "any protein identifier or FASTA (required)",
                   "taxon": "optional taxon filter", "accession": "optional pick",
                   "ph": "pH for estimated net charge (default 7.0)"},
        "returns": "unified sequence-first ProteinRecord: identity, sequence + "
                   "calculated analysis, annotations, structure availability "
                   "(experimental→predicted→sequence-only), provenance. Returns "
                   "{needs_disambiguation, candidates} when a name/gene is ambiguous.",
    },
    {
        "method": "POST", "path": "/api/protein/sequence",
        "body": {"fasta": "raw FASTA or amino-acid sequence (required)",
                 "ph": "optional pH for net charge"},
        "returns": "same ProteinRecord as GET /api/protein, resolved from the "
                   "sequence via the CRC-64 → UniParc → UniProt bridge.",
    },
    {
        "method": "GET", "path": "/api/sequence_analysis",
        "params": {"acc": "protein identifier (required)", "ph": "optional pH"},
        "returns": "Stage 3 calculated metrics: MW, pI, charge@pH, GRAVY, "
                   "composition, hydropathy profile, low-complexity regions.",
    },
    {
        "method": "GET", "path": "/api/structure_availability",
        "params": {"acc": "UniProt accession (required)"},
        "returns": "Stage 4 hierarchy: experimental PDB entities, AlphaFold "
                   "predicted models (pLDDT-gated docking suitability), or "
                   "sequence-only; with warnings. Predicted never shown as experimental.",
    },
    {
        "method": "GET", "path": "/api/homologs",
        "params": {"acc": "protein identifier (required)"},
        "returns": "structurally-characterized homologs via RCSB sequence search, "
                   "each with sequence identity and a docking caveat. Homologous "
                   "experimental structures rank above predicted models. Also "
                   "reachable via GET /api/protein?...&homologs=1.",
    },
    {
        "method": "GET", "path": "/api/domains",
        "params": {"acc": "protein identifier (required)"},
        "returns": "curated (UniProt) domains + optional InterPro families/domains "
                   "(env-gated via SNACLEX_ENABLE_INTERPRO), each source-labelled "
                   "and kept separate.",
    },
    {
        "method": "GET", "path": "/api/model_confidence",
        "params": {"acc": "UniProt accession (required)"},
        "returns": "confidence-aware AlphaFold model summary from real per-residue "
                   "pLDDT: confidence bands, low-confidence/disordered regions, and "
                   "a docking gate (low-confidence models are not marked dockable).",
    },
    {
        "method": "GET", "path": "/api/variant",
        "params": {"q": "variant, e.g. 'BRAF V600E' / 'TP53 R175H' / 'P15056:p.Val600Glu' (required)"},
        "returns": "residue-level variant analysis: sequence-coordinate mapping "
                   "with WT-residue validation, coding consequence, domain "
                   "disruption, curated-annotation overlap, known-variant match, "
                   "and clinical interpretation shown WITH evidence + review status "
                   "and limitations — never a clinical recommendation.",
    },
    {
        "method": "GET", "path": "/api/evidence",
        "params": {"protein": "protein identifier (required)",
                   "chemical": "chemical name/CID/InChIKey (required)"},
        "returns": "typed A–F protein–chemical evidence, never merged: Level B "
                   "from PubChem BioAssay (incl. inactives), plus a separate "
                   "docking (Level E) note. Docking scores are fit scores, not "
                   "affinities. Summary reports direct vs predicted counts.",
    },
    {
        "method": "GET", "path": "/api/version",
        "params": {},
        "returns": "name, version, research_only",
    },
    {
        "method": "POST", "path": "/api/jobs",
        "body": {"kind": "'dock' | 'screen' | 'benchmark'",
                 "params": "kind-specific params object"},
        "returns": "202 with {job_id, status}; poll GET /api/jobs/{id}",
        "notes": "dock: {pdb, chem, comp|pocket}. screen: {pdb, chems, comp|pocket}. "
                 "benchmark: {pdb, ligand|comp} — redocks the known ligand and "
                 "reports pocket recovery, pose RMSD, interactions recovered, and "
                 "physical plausibility.",
    },
    {
        "method": "GET", "path": "/api/benchmark/cases",
        "params": {},
        "returns": "curated known protein–ligand cases for Benchmark Mode",
    },
    {
        "method": "GET", "path": "/api/jobs/{id}",
        "params": {"id": "job id from POST /api/jobs"},
        "returns": "{status: queued|running|done|error, result?, error?}",
    },
    {
        "method": "POST", "path": "/api/upload",
        "body": "raw PDB or mmCIF text (Content-Type text/plain)",
        "returns": "upload_id + same shape as /api/analyze",
        "notes": "PDB and mmCIF formats accepted. Conservation needs a real PDB "
                 "id, so it is unavailable for uploads.",
    },
    {
        "method": "GET", "path": "/api/docs",
        "params": {},
        "returns": "this contract",
    },
]


def contract() -> dict:
    return {
        "tool": "SnaCleX",
        "version": __version__,
        "research_only": True,
        "base_url": "/",
        "limits": _LIMITS,
        "errors": "JSON {error: message} with an appropriate HTTP status "
                  "(400 bad request, 404 not found, 429 rate-limited, "
                  "500 internal, 502 upstream).",
        "endpoints": _ENDPOINTS,
    }
