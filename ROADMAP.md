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

- [x] **Add a test suite** (`tests/`, stdlib `unittest` — no pytest needed).
  Cover the pure-compute modules first since they have no network dependency:
  `pdbparse`, `pockets`, `docking`, `interactions`, `evolution`. Use the
  validation cases already named in the README (1HSG, 1CA2 → top pocket recovers
  the true ligand site) as golden tests. _(Done: 66 offline tests in `tests/`.)_
- [x] **Mock the network layer** — `snaclex/http_util.py` is the single choke
  point for all outbound fetches (RCSB, PubChem, ChEMBL, Pfam). `test_http_util`
  establishes the monkeypatch seam (retry/backoff, 4xx fast-fail, rate-limit),
  and `test_rcsb` uses it to test a client fully offline.
- [x] **GitHub Actions CI** — `.github/workflows/ci.yml` byte-compiles the
  sources and runs the suite on Python 3.11/3.12/3.13 for every push and PR.
- [x] **CHANGELOG.md** added; version now shown in the UI footer.
  _Remaining:_ a `/api/version` endpoint (deferred to Phase 2 alongside the API
  contract work).

---

## Phase 1 — Security & privacy baseline (audit priority: HIGH) — ✅ done

The audit's highest-confidence risk class. Implemented in `server.py`
(`_common_headers`, `RateLimiter`, `clean_text`) with coverage in
`tests/test_security.py`.

### 1a. Security headers (stdlib, low effort) — ✅
`_common_headers()` is emitted by both `_send_json` and `_serve_static`:
- [x] `Content-Security-Policy` — scripts locked to self + `https://3Dmol.org`
  (and `blob:` for its surface worker); images to self + PubChem; framing and
  plugins blocked. No inline scripts; inline styles allowed for `style=` attrs.
- [x] `X-Content-Type-Options: nosniff`, `Referrer-Policy: no-referrer`,
  `X-Frame-Options: DENY` (plus `frame-ancestors 'none'` in CSP).
- [x] `Strict-Transport-Security` emitted when `X-Forwarded-Proto: https`.
- [x] **CORS** deliberately same-origin (no `Access-Control-Allow-Origin`),
  documented in `_common_headers`.

### 1b. Input validation hardening (stdlib) — ✅
- [x] `clean_text()` NUL/control-char strips and length-caps `q`, `chem`, and
  `chems`; batch-screen tokens are individually bounded (`MAX_QUERY_LEN`).
- [x] XSS: added an `esc()` HTML-escaper in `app.js` for upstream-derived values,
  backed by the CSP. _(A full `innerHTML`→`textContent` sweep of all 19 sites
  remains a follow-up; CSP is the primary mitigation in the meantime.)_

### 1c. Resource-abuse / DoS controls (stdlib) — ✅
- [x] **Per-IP token-bucket rate limiter** (`RateLimiter`) on all `/api/*`, with
  a stricter budget for the expensive endpoints, returning `429` + `Retry-After`.
- [x] **Global concurrency cap** (`BoundedSemaphore`) → `503` + `Retry-After`
  when saturated. All limits tunable via env vars.
- [ ] Per-request time/size budgets for the docking search itself _(deferred to
  Phase 3, where docking moves to the async job model)_.

### 1d. Privacy & policy pages (content, low effort) — ✅
- [x] `web/privacy.html` and `web/terms.html` (no accounts/cookies/trackers,
  what's collected, retention, third-party sources + `NOTICE`, research-only
  disclaimer), linked from a new site footer.

---

## Phase 2 — Scientific provenance & reproducibility (audit priority: HIGH) — ✅ mostly done

The audit's single "highest-leverage product idea" is a *benchmark-first* mode.
**Good news: this is partly built already** — `server.py:_methods_block` emits
tool version, receptor prep, box geometry, scoring model, search seeds + random
seed, and interaction cutoffs on every dock/screen run. Extend rather than start.

- [x] **Structured export** — added alongside the existing `.txt`: **JSON**
  (full machine-readable session: metadata + every analysis + its
  methods/provenance), **CSV** (batch-screen ranking), and **PDB** (docked pose,
  from `docking.pose_to_pdb`). Buttons live in the Report tab.
  _Remaining:_ true **SDF** with bond perception (atoms-only SDF isn't useful);
  deferred to the Phase 4 docking-stack work where ligand bonds are available.
- [x] **Method cards in the UI** — `/api/pockets` and `/api/evolution` now return
  a provenance block (`snaclex/provenance.py`: method family + version + real
  parameters + scoring + interpretation + limitations), rendered as a card by
  `provenanceCardHTML`. Docking/screening already had `_methods_block`.
- [ ] **Benchmark metadata** — add a `last_benchmark_date` + dataset list to the
  methods/provenance blocks once Phase 4 benchmarks run, so each result links to
  its validation provenance.

---

## Phase 3 — Async job architecture & caching (audit priority: HIGH, effort M–L) — ✅ done

Previously work ran synchronously in the request path; caches were in-memory,
bounded LRU dicts (`_CACHE`, `_POCKET_CACHE`, `_EVO_CACHE`, `_UNIPROT_CACHE` in
`server.py`), lost on restart and not shared across instances.

**Stdlib-first (no new deps):**
- [x] **In-process job queue** (`snaclex/jobs.py`) — `/api/dock` and `/api/screen`
  now go through `POST /api/jobs` → poll `GET /api/jobs/{id}`, run on a bounded
  `ThreadPoolExecutor`, and are TTL-GC'd. Frontend uses a `submitJob()` helper
  that hides the polling. A burst queues instead of being rejected with 503.
- [x] **Disk-backed cache** (`snaclex/cache.py`) — opt-in (`SNACLEX_HTTP_CACHE`,
  on by default in `render.yaml`) TTL'd, size-bounded cache of upstream
  responses, wired into `http_util`. Stores opaque bytes only (no pickle),
  atomic writes — so restarts and repeat lookups are cheap.
- [x] **Docking grid cached per (structure, site)** (`server._get_grid`) and
  reused across every ligand in a screen and across repeat docks into the
  same site.

**Optional heavyweight track:** if `celery`/`redis` are present (env-gated),
route jobs to external workers instead of the in-process pool. The
`JobManager` submit/status interface is designed so this can be swapped in
without touching the handlers. Default deploy ignores this path.

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
