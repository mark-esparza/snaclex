# Changelog

All notable changes to SnaCleX are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### API contract & custom structures (Phase 6)
- **Documented API** — `snaclex/apidocs.py` serves a machine-readable contract at
  `GET /api/docs` (params, request bodies, limits, error shapes), rendered as a
  human page at `/api.html` and linked from the footer.
- **Custom structure upload** — `POST /api/upload` accepts a user-supplied
  PDB-format file (size-capped, validated, kept only in a bounded in-memory
  cache — never written to disk). It returns an upload id usable anywhere a PDB
  id is, so uploaded structures flow through analyze / interactions / pockets /
  docking / screening. mmCIF is detected and rejected with a clear message
  (not yet supported); conservation needs a real PDB id so it's unavailable for
  uploads. Frontend gets an "Upload PDB…" control.

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
- **Test suite** (`tests/`) — 108 stdlib `unittest` cases covering the
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
