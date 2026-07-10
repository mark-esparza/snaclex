# 07 — API Endpoint Definitions

Versioned under `/api/v1` in the target system (FastAPI + OpenAPI). The existing
unversioned endpoints (`/api/analyze`, `/api/interactions`, `/api/chemical`,
`/api/pockets`, `/api/evolution`, `/api/search`, `/api/jobs`, `/api/upload`)
remain as the **Structural Analysis** surface and are aliased under
`/api/v1/structure/*`. New platform endpoints are added; nothing is removed.

Conventions: JSON in/out; every response carries `provenance` and
`tool_version`; long work returns `202` + a `job_id`; errors are
`{error, code, source?}`. `⟳` marks endpoints that may enqueue an async job.

## Protein

| Method | Path | Purpose | Key params |
|---|---|---|---|
| GET | `/api/v1/protein/interpret` | Stage-1 query classification + candidates | `q` |
| GET | `/api/v1/protein/resolve` | Stage-2 identifier normalization + xref graph | `q`, `taxon?`, `accession?` |
| GET | `/api/v1/protein` | Full `ProteinRecord` (sequence-first) | `q` or `accession`, `isoform?` |
| POST | `/api/v1/protein/sequence` | Submit raw FASTA → resolve + analyze | body: `{fasta}` |
| GET | `/api/v1/protein/{cid}/sequence-analysis` | Stage-3 sequence metrics | `ph?` (charge at pH) |
| GET | `/api/v1/protein/{cid}/structures` | Stage-4 structure-availability hierarchy | — |
| POST | `/api/v1/protein/batch` ⟳ | Batch protein analysis | body: `{queries:[...]}` |

## Structure (existing Structural Analysis, aliased)

| Method | Path | Purpose |
|---|---|---|
| GET | `/api/v1/structure/analyze` | metadata, chains, components, raw structure (`= /api/analyze`) |
| GET | `/api/v1/structure/interactions` | atomic interaction profile (`= /api/interactions`) |
| GET | `/api/v1/structure/pockets` | LIGSITE cavities (`= /api/pockets`) |
| GET | `/api/v1/structure/evolution` | conservation (`= /api/evolution`) |
| GET | `/api/v1/structure/availability` | per-protein structure hierarchy + docking suitability |
| POST | `/api/v1/structure/upload` | analyze uploaded PDB/mmCIF (`= /api/upload`) |

## Chemical

| Method | Path | Purpose | Params |
|---|---|---|---|
| GET | `/api/v1/chemical` | PubChem properties + druglikeness (`= /api/chemical`) | `q` |
| GET | `/api/v1/chemical/normalize` | Stage-7 CID/InChIKey/parent/stereo/salt classification | `q` |
| GET | `/api/v1/chemical/{cid}/bioassays` | PubChem BioAssay evidence (Level B feed) | `cid` |
| GET | `/api/v1/chemical/search` | substructure / similarity / name search | `smiles?`, `q?`, `type` |

## Evidence & docking

| Method | Path | Purpose |
|---|---|---|
| GET | `/api/v1/evidence` | Typed A–F evidence for `protein`×`chemical` (never merged) |
| POST | `/api/v1/docking/jobs` ⟳ | Submit docking (Level E) — gated on `docking_suitable` (`= POST /api/jobs kind=dock`) |
| GET | `/api/v1/docking/jobs/{id}` | Poll docking status/result (`= /api/jobs/{id}`) |

## Variant / pathway / disease

| Method | Path | Purpose |
|---|---|---|
| GET | `/api/v1/variant` | Parse `BRAF V600E` / `acc:p.X#Y`; map to residue; consequence + evidence (no clinical call) |
| GET | `/api/v1/protein/{cid}/pathways` | Pathways (Reactome, Phase 3) |
| GET | `/api/v1/protein/{cid}/interactions` | Protein–protein interactions (IntAct/STRING, Phase 3) |

## Search & comparison

| Method | Path | Purpose |
|---|---|---|
| GET | `/api/v1/search/protein` | name/accession search |
| POST | `/api/v1/search/sequence` ⟳ | sequence-similarity search |
| GET | `/api/v1/search/ligand` | search structures by bound ligand |
| POST | `/api/v1/compare/proteins` ⟳ | family / ortholog / WT-vs-mutant comparison |
| POST | `/api/v1/compare/targets` ⟳ | one compound across many proteins (and vice-versa) |

## Report / provenance / export

| Method | Path | Purpose |
|---|---|---|
| GET | `/api/v1/report/{cid}` | Assembled 13-tab report payload |
| POST | `/api/v1/report/{cid}/export` | Downloadable provenance-aware snapshot (JSON/CSV/PDF) |
| GET | `/api/v1/provenance/{evidence_id}` | Inspect a single evidence object |
| GET | `/api/v1/version` · `/api/v1/docs` | version · OpenAPI (`= /api/version`, `/api/docs`) |

## Async job envelope (shared)

`POST …/jobs → 202 {job_id}`; `GET …/jobs/{id} → {status: queued|running|done|error, result?, error?}`.
Used for docking, large alignments, batch, comparison, and (optional)
prediction — never inline (NFR-4). Backed today by `snaclex/jobs.py`
(ThreadPoolExecutor), swappable for Celery/Redis.

## Vertical-slice endpoints shipping now (stdlib `server.py`)

Unversioned to match the current app; they are the `/api/v1` shapes above minus
the version prefix:

- `GET /api/protein?q=…[&taxon=…][&accession=…]` — interpret → resolve → build
  `ProteinRecord` (UniProt/UniParc/NCBI/RCSB/AlphaFold) + sequence analysis +
  structure availability. Degrades to sequence-only.
- `GET /api/resolve?q=…` — Stage-1/2 only (classification + xref graph).
- `GET /api/structure_availability?acc=…` — Stage-4 hierarchy for a UniProt acc.
- `GET /api/evidence?protein=…&chemical=…` — typed A–F evidence (Level A from PDB
  bound ligands, Level B from PubChem BioAssay; C/D/E/F scaffolded with clear
  "not yet gathered" markers where a source is Phase 2/3).
- `GET /api/sequence_analysis?acc=…[&ph=7.4]` or `POST /api/protein/sequence`
  with `{fasta}` — Stage-3 pure-compute metrics.

All are added to `apidocs.py` so `GET /api/docs` stays the single source of
truth, and all reuse the existing security/rate-limit/observability pipeline in
`server.py`.
