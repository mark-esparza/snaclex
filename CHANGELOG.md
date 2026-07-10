# Changelog

All notable changes to SnaCleX are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Sequence-first platform (Phase 1 — expansion)
- **New capability: analyze any cataloged protein, structure or not.** SnaCleX
  now builds a unified, sequence-first **ProteinRecord** that exists even when no
  crystal/cryo-EM/NMR structure is available — closing the app's biggest gap
  (previously PDB-only). Full design set in [`docs/platform/`](docs/platform/).
- **New source adapters** (separate, offline-testable, provenance-stamped):
  `uniprot.py` (UniProtKB — primary annotation, Swiss-Prot vs TrEMBL),
  `uniparc.py` (UniParc — sequence archive / identity bridge),
  `ncbi.py` (NCBI Protein / RefSeq via E-utilities, optional `NCBI_API_KEY`),
  `alphafold.py` (AlphaFold DB — predicted models, pLDDT bands, PAE, docking gate).
  RCSB and PubChem adapters extended (by-UniProt structure search; BioAssay
  evidence + InChIKey/parent normalization).
- **New services**: `idresolve.py` (query interpretation + cross-reference graph;
  identity by accession+taxon+checksum, never by name/gene alone),
  `seqanalysis.py` (calculated MW/pI/charge@pH/GRAVY/composition/low-complexity —
  pure, always available), `proteinrecord.py` (record assembly + Stage-4 structure
  hierarchy), `evidence.py` (**universal evidence object + A–F levels with hard
  invariants** — a docking score can never be an affinity; levels never merge),
  `checksum.py` (SWISS-PROT CRC-64 + MD5/SHA-256 sequence identity).
- **New API endpoints** (documented in `/api/docs`): `GET /api/resolve`,
  `GET /api/protein`, `POST /api/protein/sequence`, `GET /api/sequence_analysis`,
  `GET /api/structure_availability`, `GET /api/evidence`.
- **Existing structural tooling preserved** as the modular Structural Analysis
  component (analyze/interactions/pockets/evolution/docking unchanged).
- **~60 new offline tests** (adapters, evidence invariants, sequence analysis,
  identifier resolution, record assembly, platform endpoints). Suite stays fully
  offline and dependency-free. Live source fields are marked `⚠ VERIFY` where the
  environment's egress policy blocked round-trip confirmation.

### Predicted structures & enrichment (Phase 2 — in progress)
- **Confidence-aware AlphaFold handling** — `alphafold.build_confidence` derives
  confidence bands, low-confidence and disordered regions, and a docking gate
  from the model's **real per-residue pLDDT** (B-factor column), not just the API
  summary metric. New `GET /api/model_confidence?acc=...` fetches coordinates on
  demand (kept out of the lightweight availability path) so a low-confidence
  model is never marked dockable.
- **InterPro enrichment adapter** (`interpro.py`) — optional and env-gated
  (`SNACLEX_ENABLE_INTERPRO`), provenance-stamped, degrades gracefully; families
  and domains are kept **separate** from curated UniProt domains, each
  source-labelled. New `GET /api/domains?acc=...`.
- +10 offline tests (InterPro parse/gating, confidence gate, both endpoints).
- **Homology-based structure selection** — `rcsb.sequence_search` (RCSB mmseqs2)
  + `homology.py` find structurally-characterized homologs when a protein has no
  direct structure. The Stage-4 hierarchy now ranks **experimental → homologous
  experimental → predicted → sequence-only**; a homolog is only offered as a
  docking receptor above a 0.40 identity floor and always carries a "not the
  requested protein" caveat. New `GET /api/homologs?acc=...` and
  `GET /api/protein?...&homologs=1`. The evidence endpoint surfaces Level-D
  homolog *candidates* (`?homologs=1`) without asserting a "binds" claim. +13
  offline tests.

### Observability (Phase 7b)
- **Structured request logging** — every request logs one line to stdout
  (`METHOD path -> status durationms ip=<hash>`) via the stdlib `logging` module,
  so the public service finally has access logs (Render captures stdout). Client
  IPs are logged only as a **salted hash** (per the privacy policy); job-status
  polls are demoted to DEBUG. `SNACLEX_LOG_LEVEL` / `SNACLEX_IP_SALT` tune it.
- **5xx tracebacks logged server-side**; the client now gets a generic
  "Internal error — please retry." instead of the raw exception text (closes the
  audit's info-leak finding). Failed jobs are logged from the worker too.

### Changed
- **Product positioning** — reframed the tagline, page title/description, and
  README as *"a reproducible, browser-based structural-biology workbench"* for
  interaction analysis, pocket prioritization, conservation-guided
  interpretation, and exploratory docking. Explicitly positioned as an
  interpretation workbench, **not** a full drug-discovery platform, with a
  "Where it fits" comparison (vs. RCSB viewing, command-line docking, basic
  viewers, and black-box AI docking demos).

### Testing & CI
- **Browser end-to-end smoke test** (`e2e/`, Playwright) — boots the real server
  and drives Chromium through load → pick a bound molecule → interactions, and
  **asserts no uncaught JS error** fires (the audit's P0: the front end had only
  ever been syntax-checked). Runs fully offline via a pre-seeded disk cache (no
  upstream calls), in a separate CI workflow (`.github/workflows/e2e.yml`) so the
  core stdlib suite stays dependency-free. Playwright is a dev/CI-only dep
  (`requirements-dev.txt`). Plus a pre-deploy **manual QA checklist**
  (`docs/qa-checklist.md`).

### Added
- **Benchmark Mode** — a new tab and `benchmark` job kind that runs SnaCleX
  against known protein–ligand cases and reports research-credibility metrics:
  **pocket recovered** (yes/no + distance/rank vs the crystal site), **pose
  RMSD** (re-docking the known ligand into its own site), **interactions
  recovered** (Y/Z of the experimental contact residues), and **physical
  plausibility** (pass/fail clash check). Curated cases (1HSG/indinavir,
  3PTB/benzamidine, 4DFR/methotrexate) are served from `GET /api/benchmark/cases`;
  any loaded structure's ligand can also be benchmarked. (`benchmark.benchmark_case`)
- **Per-chain loading for large assemblies** — when a structure exceeds the
  interactive atom limit, `/api/analyze` returns its chain list
  (`{too_large, n_atoms, limit, chains:[{chain, atom_count}]}`) and the UI shows
  a chain picker. `/api/analyze?pdb=ID&chain=X` loads just that chain (subset
  cached under a synthetic id, usable across all analyses). This is how big
  mega-assemblies like `6BCX` (113,610 atoms) become analyzable one chain at a
  time. (`pdbparse.subset_chain`)

### Fixed
- **Structures with no legacy PDB file now load** (e.g. `6BCX`, which 404'd on
  `files.rcsb.org/download/*.pdb`). `rcsb.fetch_structure` falls back to the
  mmCIF (`.cif`) file, and a new `pdbparse.parse_mmcif` / `parse_structure`
  auto-detects and parses the `_atom_site` loop. mmCIF uploads are accepted too.
- **mmCIF structures render in the viewer** — the 3Dmol viewer only reads PDB
  format, so the server now re-serializes mmCIF-sourced structures to PDB
  (`pdbparse.to_pdb`) before sending them. Previously they were handed to 3Dmol
  as `"pdb"` and failed to display.
- **Oversized structures no longer freeze the page** — entries above the
  PDB-format / interactive limit (99,999 atoms; `SNACLEX_MAX_ATOMS`) are
  rejected with a clear `413` message instead of timing out the tab. mmCIF-only
  mega-assemblies are exactly this case.

### API contract & custom structures (Phase 6)
- **Documented API** — `snaclex/apidocs.py` serves a machine-readable contract at
  `GET /api/docs` (params, request bodies, limits, error shapes), rendered as a
  human page at `/api.html` and linked from the footer.
- **Custom structure upload** — `POST /api/upload` accepts a user-supplied
  PDB-format file (size-capped, validated, kept only in a bounded in-memory
  cache — never written to disk). It returns an upload id usable anywhere a PDB
  id is, so uploaded structures flow through analyze / interactions / pockets /
  docking / screening. Both PDB and mmCIF are accepted; conservation needs a
  real PDB id so it's unavailable for uploads. Frontend gets an "Upload PDB…"
  control.

### Accessibility & onboarding (Phase 5)
- **Example structures** — one-click 1HSG / 1CA2 / 4HHB buttons in the load card.
- **Accessibility** — visible keyboard-focus styles (`:focus-visible`), an
  `aria-live` status region, and labelled inputs.
- **Responsive** — the single-screen layout stacks on narrow viewports
  (≤ 820 px) instead of crowding.

### Docking rigor (Phase 4, partial)
- **Redocking benchmark harness** (`snaclex/benchmark.py` + `python -m
  snaclex.benchmark`) — self-docks each crystallographic ligand back into its
  receptor and reports pose-recovery metrics (top-1 RMSD ≤ 2 Å, median/mean).
  Dependency-free, runs on local PDB files or fetched IDs; scales to
  PoseBusters/CrossDocked/PDBbind. This is the audit's "measure the docker
  before changing it" step.
- **Benchmark metadata in the methods block** (the deferred Phase 2 item) — if
  `benchmark_results.json` is committed, every docking result links to the
  method's last measured numbers (`provenance.docking_benchmark`).
- _Deferred:_ the optional AutoDock Vina / GNINA upgrade track (needs those
  engines + RDKit/Meeko installed) — out of scope for the zero-dependency
  default; the same `summarize` will benchmark it head-to-head once present.

### Performance & scaling (Phase 3)
- **Async job queue** (`snaclex/jobs.py`) — docking and batch screening now run
  on a bounded `ThreadPoolExecutor` via `POST /api/jobs` → poll
  `GET /api/jobs/{id}`, instead of blocking the HTTP request. A burst now
  *queues* behind the worker pool rather than being rejected, and slow computes
  no longer risk proxy/browser timeouts. Jobs are TTL'd and garbage-collected.
- **Disk-backed HTTP cache** (`snaclex/cache.py`) — opt-in (`SNACLEX_HTTP_CACHE`)
  TTL'd, size-bounded cache of upstream (RCSB/PubChem/ChEMBL/Pfam) responses,
  so repeat lookups are instant and survive restarts. Stores opaque bytes only
  (no pickling), writes atomically. Enabled by default in `render.yaml`.
- **Docking-grid cache** — the scoring grid is cached per (structure, site) and
  reused across every ligand in a screen and across repeat docks into the
  same site.
- The frontend submits these via a `submitJob()` helper that hides the polling,
  so the render code is unchanged.

### Provenance & reproducibility (Phase 2)
- **Method-transparency blocks** for pocket detection and conservation
  (`snaclex/provenance.py`) — method family, version, real parameters (pulled
  from the module constants), scoring formula, a plain-language "what this score
  means / does not mean", and limitations. Now returned by `/api/pockets` and
  `/api/evolution`, matching the methods block docking/screening already emit.
- **Method/interpretation cards** rendered in the Pockets and Evolution tabs
  (`provenanceCardHTML` in `app.js`).
- **Structured exports** beyond the existing `.txt`:
  - **JSON** — the full machine-readable session (metadata + every analysis +
    its methods/provenance) for reproducibility and audit.
  - **CSV** — the batch-screen ranking, spreadsheet/pandas-ready.
  - **PDB** — the docked pose as a real coordinate file for PyMOL/Chimera/etc.
- **`/api/version`** endpoint (the item deferred from Phase 0).

### Security (Phase 1)
- **Security headers** on every response — a tuned `Content-Security-Policy`
  (locks scripts to self + the 3Dmol.js CDN, images to self + PubChem, blocks
  framing and plugins), plus `X-Content-Type-Options`, `Referrer-Policy`,
  `X-Frame-Options`, and `Strict-Transport-Security` when served over HTTPS.
  CORS is deliberately same-origin (no `Access-Control-Allow-Origin`).
- **Abuse / DoS controls** — a thread-safe per-IP token-bucket rate limiter on
  all `/api/*` calls, a stricter budget for the compute-heavy endpoints
  (`/dock`, `/screen`, `/pockets`, `/evolution`), and a global concurrency cap
  that returns `503` with `Retry-After` when saturated. All tunable via env vars.
- **Input hardening** — free-text query params are NUL/control-char stripped and
  length-capped; batch-screen tokens are individually bounded.
- **XSS defense-in-depth** — upstream-derived values interpolated into the DOM
  are HTML-escaped (`escapeHtml`), complementing the CSP.

### Added
- **Privacy & Terms pages** (`web/privacy.html`, `web/terms.html`) linked from a
  new site footer: no accounts/cookies/trackers, what's collected, third-party
  data sources, retention, and the research-only disclaimer.
- **Test suite** (`tests/`) — 128 stdlib `unittest` cases covering the
  pure-compute core (`pdbparse`, `interactions`, `docking`, `pockets`,
  `report`) plus the PubChem/RCSB helpers. Runs fully offline; the HTTP fetch
  layer (`http_util`) is exercised via a mock seam (retry/backoff, fast-fail on
  4xx, rate-limit handling), establishing the pattern for testing the network
  clients without live API calls.
- **Continuous integration** (`.github/workflows/ci.yml`) — byte-compiles all
  sources and runs the test suite on Python 3.11/3.12/3.13 for every push and
  pull request.
- **ROADMAP.md** — prioritized plan derived from the external audit report,
  mapped to the actual codebase and the project's zero-dependency philosophy.

_This covers Phase 0 (engineering hygiene), Phase 1 (security & privacy
baseline), Phase 2 (provenance & structured export), Phase 3 (async jobs & caching),
Phase 4 (docking benchmark harness), Phase 5 (accessibility & onboarding), and
Phase 6 (API contract & structure upload) of the roadmap._

## [0.1.0]

- Initial public research tool: PDB structure loading (RCSB), atomic
  interaction profiling, LIGSITE pocket detection, Pfam-based conservation
  scoring, pure-Python Monte-Carlo docking, batch screening, PubChem/ChEMBL
  chemical lookup, and `.txt` session report export.
