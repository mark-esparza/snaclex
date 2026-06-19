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
        "method": "GET", "path": "/api/interface",
        "params": {"pdb": "structure id (required)",
                   "a": "comma-separated chain ids for group A (optional)",
                   "b": "comma-separated chain ids for group B (optional)"},
        "returns": "protein-protein/peptide interface profile (epitope/paratope "
                   "contact residues) + methods/provenance. Defaults to the two "
                   "smallest chains vs the rest when a/b are omitted.",
    },
    {
        "method": "GET", "path": "/api/nucleic",
        "params": {"pdb": "structure id (required)",
                   "chain": "optional nucleic chain to restrict to"},
        "returns": "protein-DNA/RNA interface profile (binding residues + "
                   "contacted nucleotides) + methods/provenance",
    },
    {
        "method": "GET", "path": "/api/hla",
        "params": {"pdb": "structure id (required)"},
        "returns": "HLA/MHC detection, peptide-binding-groove + anchor-pocket "
                   "residues, peptide contacts + methods/provenance",
    },
    {
        "method": "GET", "path": "/api/hla/cases",
        "params": {},
        "returns": "curated well-known HLA structures (allele, pdb, class, note)",
    },
    {
        "method": "GET", "path": "/api/esm",
        "params": {},
        "returns": "ESM model-layer config/status (available, fold/score models). "
                   "Features are active only when ESM_API_KEY is set.",
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
        "method": "GET", "path": "/api/version",
        "params": {},
        "returns": "name, version, research_only",
    },
    {
        "method": "POST", "path": "/api/jobs",
        "body": {"kind": "'dock' | 'screen' | 'benchmark' | 'variants'",
                 "params": "kind-specific params object"},
        "returns": "202 with {job_id, status}; poll GET /api/jobs/{id}",
        "notes": "dock: {pdb, chem, comp|pocket}. screen: {pdb, chems, comp|pocket}. "
                 "benchmark: {pdb, ligand|comp} — redocks the known ligand and "
                 "reports pocket recovery, pose RMSD, interactions recovered, and "
                 "physical plausibility. variants: {pdb, variants:[...], uniprot?} — "
                 "maps protein-position variants onto residues (pocket / "
                 "conservation / HLA-groove context), TOPMed/BRAVO frequency for "
                 "genomic inputs, and optional ESM scores.",
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
