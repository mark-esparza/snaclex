# 08 — Backend Repository Structure

Two layouts: the **current + slice** layout (what exists after this change) and
the **target** layout (the FastAPI/Postgres destination). The slice extends the
flat `snaclex/` package so imports and the stdlib server keep working; the target
reorganizes into a `src/` service tree once the heavyweight stack is adopted.

## 8.1 Current + slice (ships now)

```
snaclex/
  __init__.py            version
  http_util.py    ✓      shared HTTP: retry/backoff/cache/rate-limit choke point
  cache.py        ✓      TTL disk cache
  jobs.py         ✓      async job queue (ThreadPool; Celery-swappable)
  provenance.py   ✓      method-transparency blocks  (+ new: record/evidence provenance)
  apidocs.py      ✓      machine-readable API contract (+ new endpoints)

  # ── Source adapters (one per source) ───────────────────────────
  rcsb.py         ✓      RCSB PDB
  pubchem.py      ✓      PubChem (+ new: bioassay_summary, normalize/InChIKey/parent)
  chembl.py       ✓      optional drug/bioactivity cross-ref
  uniprot.py      🆕     UniProtKB
  uniparc.py      🆕     UniParc
  ncbi.py         🆕     NCBI Protein / RefSeq (E-utilities)
  alphafold.py    🆕     AlphaFold DB

  # ── Structural Analysis component (existing, preserved) ────────
  pdbparse.py     ✓      PDB + mmCIF parser
  interactions.py ✓      atomic interaction profiler
  pockets.py      ✓      LIGSITE cavity finder
  docking.py      ✓      grid-map Monte-Carlo docker
  evolution.py    ✓      Pfam conservation
  benchmark.py    ✓      redocking benchmark CLI

  # ── Platform services (new, sequence-first) ───────────────────
  seqanalysis.py  🆕     Stage-3 pure-compute sequence metrics
  idresolve.py    🆕     Stage-1/2 query interpretation + xref graph
  proteinrecord.py 🆕    unified ProteinRecord assembly + structure availability
  evidence.py     🆕     universal evidence object + A–F levels + invariants

server.py         ✓      stdlib server; new routes wired in
web/              ✓      vanilla-JS SPA (+ new protein-report views, Phase 1.5)
tests/            ✓      offline unittest (+ new: seqanalysis, evidence, adapters, idresolve)
docs/platform/    🆕     this design set
```

Design rule: adapters return normalized fragments + provenance and never leak raw
source JSON upward; services (`idresolve`, `proteinrecord`, `evidence`,
`seqanalysis`) compose adapter outputs and contain no HTTP. This keeps services
unit-testable offline (feed them recorded fragments) exactly like the existing
compute modules.

## 8.2 Target layout (FastAPI / Postgres destination)

```
backend/
  pyproject.toml                # deps: fastapi, pydantic, sqlalchemy, alembic, httpx, redis, celery, boto3
  src/snaclex/
    api/                        # FastAPI routers, one per resource (protein, chemical, evidence, ...)
      v1/ __init__.py protein.py structure.py chemical.py evidence.py variant.py search.py report.py
      deps.py  errors.py  openapi.py
    core/                       # config, logging, rate-limit, api-key mgmt, provenance envelope
    adapters/                   # uniprot.py uniparc.py ncbi.py rcsb.py alphafold.py pubchem.py + optional/
      base.py                   # Adapter protocol: fetch → normalized fragment + provenance
      http.py                   # httpx client: retry, backoff, conditional GET, per-source limiter
    services/
      resolution.py sequence_analysis.py structure_analysis.py
      chemical_normalization.py evidence_integration.py search.py comparison.py
    structural/                 # the migrated existing compute modules (pdbparse, interactions, ...)
    models/                     # Pydantic (API) + SQLAlchemy (DB) + ProteinRecord/Evidence
    db/                         # session, migrations (alembic), repositories
    jobs/                       # celery app + tasks (docking, alignment, batch, prediction)
    storage/                    # object-store client (structures, alignments, report snapshots)
    graph/                      # optional graph-DB projection
  tests/  unit/  integration/  fixtures/     # pytest
  Dockerfile  docker-compose.yml
frontend/                       # React + TS + Mol* (see 09-frontend-ia.md)
clients/python/                 # thin Python client + CLI over the versioned API
```

## 8.3 Adapter interface (both layouts)

```python
# base contract every source adapter honors
def fetch(identifier: str, *, ctx) -> Fragment:
    """Return a normalized fragment + provenance, or a `unavailable` fragment.
    Never raises for a routine miss; never returns raw upstream JSON."""

Fragment = {
    "data": {...},                 # normalized, source-agnostic keys
    "provenance": {                # the envelope from 04-data-model §4.0
        "source": "UniProtKB", "source_id": "...", "source_version": "...",
        "retrieved_utc": "...", "evidence_category": "curated", "limitations": None,
    },
    "status": "ok" | "unavailable" | "rate_limited",
}
```

Composition services import adapters, not `http_util`; only adapters do I/O. That
boundary is what makes the whole platform testable offline and lets any single
source be disabled without touching callers (NFR-2).
