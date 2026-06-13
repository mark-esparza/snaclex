# SnaCleX Roadmap

This roadmap turns the findings of the external **SnaCleX audit report** into a
concrete, prioritized plan of work, mapped against the *actual* codebase (the
audit was written from the public UI only, with no source access).

Each item notes what the audit recommended, what the code **already does**, the
**real gap**, and concrete next steps with file-level pointers.

## Guiding constraint: the zero-dependency philosophy

SnaCleX is intentionally **pure Python standard library** (`requirements.txt` is
empty by design; `render.yaml` has no build step). The audit recommends a heavy
stack — FastAPI, Celery/Redis, AutoDock Vina, GNINA, RDKit, OpenTelemetry,
Sentry. Those are good *destinations*, but adopting them wholesale would discard
the project's main operational virtue: it deploys anywhere with no build and no
dependency supply chain.

This roadmap therefore splits each recommendation into:

- **Stdlib-first** — what we can do now without new dependencies (most of the
  security, provenance, UX, and export work falls here).
- **Optional heavyweight track** — opt-in, env-gated integrations (Vina/GNINA
  docking, async workers) that only activate when the extra packages are
  present, so the default deploy stays dependency-free.

---

## Phase 0 — Engineering hygiene (foundation, do first)

The audit flags "low observable production maturity": no tests, no CI, no
version/changelog surfaced. None of the later phases are safe to land without a
test harness.

- [ ] **Add a test suite** (`tests/`, stdlib `unittest` — no pytest needed).
  Cover the pure-compute modules first since they have no network dependency:
  `pdbparse`, `pockets`, `docking`, `interactions`, `evolution`. Use the
  validation cases already named in the README (1HSG, 1CA2 → top pocket recovers
  the true ligand site) as golden tests.
- [ ] **Mock the network layer** — `snaclex/http_util.py` is the single choke
  point for all outbound fetches (RCSB, PubChem, ChEMBL, Pfam). Inject a
  fixture/monkeypatch seam there so `rcsb`/`pubchem`/`chembl`/`evolution` are
  testable offline with checked-in sample payloads.
- [ ] **GitHub Actions CI** — run the test suite + `python -m py_compile` on
  push/PR. Cheap, and directly answers the audit's "no operational maturity"
  point.
- [ ] **CHANGELOG.md + version surfacing** — `__version__` exists in
  `snaclex/__init__.py` (`0.1.0`) and is already echoed in the methods block; also
  surface it in the UI footer and add a `/api/version` endpoint.

---

## Phase 1 — Security & privacy baseline (audit priority: HIGH)

The audit's highest-confidence risk class. The current `Handler` sends no
security headers and has no rate limiting; user inputs flow into a rich UI.

### 1a. Security headers (stdlib, low effort)
Add a single `_send_headers` helper used by both `_send_json` and
`_serve_static` in `server.py`:
- [ ] `Content-Security-Policy` — the UI loads 3Dmol.js; lock script sources to
  self + the specific CDN, disallow inline where feasible.
- [ ] `X-Content-Type-Options: nosniff`, `Referrer-Policy: no-referrer`,
  `X-Frame-Options: DENY` (or `frame-ancestors` in CSP).
- [ ] `Strict-Transport-Security` when served over HTTPS (Render terminates TLS).
- [ ] Explicit, narrow **CORS** policy (default: same-origin only). Today there
  is no `Access-Control-Allow-Origin` header at all — make the absence
  deliberate and documented.

### 1b. Input validation hardening (stdlib)
Current validation is partial: `rcsb.normalize_pdb_id`, integer parsing of
`comp`/`pocket`, and a 10-item cap on screening. Tighten:
- [ ] Strict allowlist/length bounds on `pdb` (4-char alphanumeric), `chem`/`q`
  (length-capped, control-char stripped), and `chems` tokens.
- [ ] Confirm all server-derived text is safely encoded where `web/app.js`
  injects it into the DOM (audit's XSS concern) — prefer `textContent` over
  `innerHTML` for any value originating from RCSB/PubChem/ChEMBL.

### 1c. Resource-abuse / DoS controls (stdlib)
The audit's most credible risk: pocketing, conservation, docking, and screening
are compute-heavy on a public anonymous service running a `ThreadingHTTPServer`.
- [ ] **Per-IP token-bucket rate limiter** in `Handler` (in-memory, stdlib) on
  the expensive endpoints (`/api/dock`, `/api/screen`, `/api/pockets`,
  `/api/evolution`), returning `429` with `Retry-After`.
- [ ] **Global concurrency cap** (bounded worker semaphore) so a burst of docking
  jobs can't exhaust the box; return `503` with backpressure when saturated.
- [ ] Per-request time/size budgets for the docking search.

### 1d. Privacy & policy pages (content, low effort)
The audit's clearest external-compliance gap: no privacy/terms/cookie notice.
- [ ] Add `web/privacy.html` and `web/terms.html`: what's collected (request
  logs/IP only — no accounts, no PII inputs), retention, third-party data
  sources (RCSB/PubChem/ChEMBL/Pfam) and their terms (already in `NOTICE`),
  and the research-only disclaimer. Link from the UI footer.

---

## Phase 2 — Scientific provenance & reproducibility (audit priority: HIGH)

The audit's single "highest-leverage product idea" is a *benchmark-first* mode.
**Good news: this is partly built already** — `server.py:_methods_block` emits
tool version, receptor prep, box geometry, scoring model, search seeds + random
seed, and interaction cutoffs on every dock/screen run. Extend rather than start.

- [ ] **Structured export** — today the report is `.txt` only
  (`snaclex/report.py`). Add **JSON** export of the full session (metadata +
  methods block + profiles + docking/screen results) and **CSV** for screen
  rankings. This is the audit's "no structured export (JSON/CSV/SDF)" gap. Add
  **SDF/PDB** export for docked poses (`docking.pose_to_pdb` already produces
  PDB text).
- [ ] **Method cards in the UI** — per analysis tab (Pockets/Evolution/Docking/
  Screening), surface a compact panel: method family + version, known
  limitations, and a "What this score means / does *not* mean" note. The
  scoring/disclaimer strings already exist in `_methods_block`; render them.
- [ ] **Benchmark metadata** — add a `last_benchmark_date` + dataset list to the
  methods block once Phase 4 benchmarks run, so each result links to its
  validation provenance.

---

## Phase 3 — Async job architecture & caching (audit priority: HIGH, effort M–L)

Today work runs synchronously in the request path; caches are in-memory, bounded
LRU dicts (`_CACHE`, `_POCKET_CACHE`, `_EVO_CACHE`, `_UNIPROT_CACHE` in
`server.py`), lost on restart and not shared across instances.

**Stdlib-first (no new deps):**
- [ ] **In-process job queue** — convert `/api/dock` and `/api/screen` to a
  submit→poll pattern: `POST` returns a `job_id`; a stdlib
  `ThreadPoolExecutor`-backed worker pool runs the job; `GET /api/jobs/{id}`
  returns status/result. Reuses the concurrency cap from Phase 1c.
- [ ] **Disk-backed cache** for upstream fetches and derived features (parsed
  structures, pockets, conservation) so restarts and cold tabs are cheap — a
  simple keyed JSON/pickle cache directory, TTL'd, stdlib only.
- [ ] **Precompute the docking grid once per (structure, site)** and reuse across
  screening ligands — `docking.build_grid` already supports this; cache the grid.

**Optional heavyweight track:** if `celery`/`redis` are present (env-gated),
route jobs to external workers instead of the in-process pool. Default deploy
ignores this path.

---

## Phase 4 — Docking rigor & benchmark suite (audit priority: MEDIUM, very high impact)

The audit's clearest *scientific* gap: the current docker is an honest,
self-disclosed "approximate, relative-score" Monte-Carlo rigid-body search
(`snaclex/docking.py`; README explains Vina/RDKit were avoided for clean
zero-dependency Windows installs).

- [ ] **Benchmark harness** (can be stdlib + checked-in test structures): measure
  the *current* docker before changing it — Top-1 RMSD ≤ 2 Å and median RMSD on
  a small redock set. `docking.rmsd_to_reference` already exists and `/api/dock`
  already computes redock RMSD when the same ligand is crystallized. Formalize
  this into a scriptable benchmark with a results table.
  - Datasets from the audit: PoseBusters (pose validity), CrossDocked2020
    (redock/cross-dock), PDBbind/BindingDB (affinity-linked), LigASite (pocket
    sites).
- [ ] **Virtual-screening metrics** — EF1%, BEDROC, AUROC on a curated
  target/decoy set for `/api/screen`.
- [ ] **Optional Vina/GNINA track** — env-gated: if AutoDock Vina (and optionally
  GNINA reranking) + Meeko prep are installed, expose them as an *upgraded*
  docking mode selectable in the UI, benchmarked head-to-head against the
  built-in docker. The pure-Python docker remains the dependency-free default.

---

## Phase 5 — Accessibility & onboarding (audit priority: MEDIUM, effort S–M)

Audit: onboarding structure is good, but no help/glossary/examples, and a11y is
unverified (color reliance, focus, ARIA, keyboard).

- [ ] **Sample/"known-good" walkthroughs** — one-click load of 1HSG / 1CA2 (the
  README validation cases) as guided examples.
- [ ] **WCAG pass on `web/`** — visible focus states, ARIA labels on viewer
  controls, non-color cues for interaction types (the viewer currently leans on
  color themes + interaction-line color), contrast check on `style.css`.
- [ ] **Interpretation overlays / glossary** — short "what this means" tooltips
  for pocket druggability score, conservation score, and relative docking score.
- [ ] **Responsive/narrow-viewport mode** — collapsible analysis panels; the
  single-screen layout is desktop-first today.

---

## Phase 6 — API contract & custom structure upload (audit priority: MEDIUM→LOWER)

- [ ] **Documented API contract** — the JSON endpoints in `server.py` are
  undocumented. Publish a reference (handwritten OpenAPI/JSON doc, served at
  `/api/docs`) covering params, limits, and error shapes. Enables scripting and
  peer audit.
- [ ] **Custom structure upload** (audit's highest-value, highest-risk item) —
  accept user-uploaded mmCIF/PDB and AlphaFold/UniProt import. `snaclex/pdbparse.py`
  already parses PDB text; the work is a validated upload path + size limits +
  the security scope it expands (untrusted file parsing, storage). Gate behind
  Phase 1 hardening.

---

## Priority summary

| Phase | Theme | Audit priority | Effort | New deps? |
|-------|-------|----------------|--------|-----------|
| 0 | Tests + CI + changelog | (foundation) | S–M | none |
| 1 | Security & privacy baseline | HIGH | M | none |
| 2 | Provenance & structured export | HIGH | M | none |
| 3 | Async jobs + durable caching | HIGH | M–L | none (opt: celery/redis) |
| 4 | Docking benchmarks + Vina/GNINA | MEDIUM | L | opt: vina/gnina/rdkit |
| 5 | Accessibility & onboarding | MEDIUM | S–M | none |
| 6 | API docs + structure upload | MEDIUM→LOWER | M–L | none |

**Recommended sequence:** Phase 0 → 1 → 2 in order (each is low-risk, high-value,
and dependency-free), then 3, then pick up 4/5/6 in parallel. Phases 0–2 and 5
keep SnaCleX entirely within its zero-dependency philosophy; only the *optional*
tracks in 3 and 4 introduce packages, and only when explicitly enabled.

---

*Derived from the external SnaCleX audit report and a review of the current
source tree (`server.py`, `snaclex/`, `web/`). Items marked "already" reflect
behavior the audit could not observe from the public UI but which exists in code.*
