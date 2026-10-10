# Changelog

All notable changes to SnaCleX are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added: upstream connectivity self-test
- **`snaclex/selftest.py` + `GET /api/selftest` + `python -m snaclex.selftest`** —
  SnaCleX has no database of its own; every result is assembled live from five
  public services, so a deployment is only as good as its outbound access.
  This checks all of them in one pass and reports per-source status, latency,
  and what each check proves.
- **The checks exercise the real integration, not just reachability.** Each one
  calls SnaCleX's own client for that source and asserts a known fact about the
  result — aspirin is CID 2244, crambin parses to a few hundred protein atoms,
  1CRN maps to a UniProt accession — so a silently changed response shape fails
  here rather than in a user's analysis. Eight checks across RCSB
  (files/data/search/GraphQL), PubChem (properties + 3D conformer), ChEMBL and
  InterPro/Pfam.
- **Exit status is 0 only when every check passes**, so it can gate a deploy.
- The endpoint is rate-limited as an expensive route and cached for 60 s, so it
  cannot be used to amplify traffic at RCSB/PubChem/EBI on our behalf.
- `http_util.fail_fast()` — a **thread-local** retry/timeout override used by the
  self-test, so a diagnostic reports a down host in ~0.3 s instead of spending
  the full retry budget, **without** altering concurrent real requests on the
  threaded server. The full sweep went from 76 s to under 4 s.

### Fixed
- `HEAD` requests returned **501** from every URL (`BaseHTTPRequestHandler`
  answers "Unsupported method" unless `do_HEAD` exists), so uptime monitors,
  load balancers and link checkers all got an error. HEAD now returns the same
  status, `Content-Length` and security headers as GET, with no body. Render's
  own health check uses GET, so this was latent rather than a live outage.


### Interface
- **Institutional seal** replaces the placeholder glyph in the masthead — an
  inline SVG mark (ring, hexagon, gold core, three nodes), so it adds no
  network origin and needs nothing from the CSP's `img-src`.
- **The section nav's dead space now carries context** — the loaded structure's
  id, protein atom count and bound-component count, right-aligned. The
  breadcrumb scrolls away on a long results page; this does not.

### Fixed
- **`HEAD` requests returned 501.** `BaseHTTPRequestHandler` answers HEAD with
  "Unsupported method" unless `do_HEAD` exists, so uptime monitors, load
  balancers and link checkers all got an error from every URL. HEAD now
  returns the same status, `Content-Length` and security headers as GET with
  no body. Render's own health check uses GET, so this was latent rather than
  a live outage.


### Interface: from desktop window to research site
The theme kept its palette, fonts and widget treatment, but the page was built
as an application *window* — a caption bar and menu bar wrapped around a tab
control, each layer with its own frame. That read as a tab inside a tab. The
architecture is now the one a public research site of the period actually used:

- **Masthead band → site nav strip → breadcrumb → sidebar + content → agency
  footer**, each band full-bleed with its contents on one shared measure. One
  chrome layer instead of four nested frames.
- **The window caption is now a masthead** (no caption buttons, no rounded
  corners, no floating-on-a-desktop shadow), and the decorative `File / Edit /
  View` menu bar is replaced by a **real site nav** — Workbench, API Reference,
  RCSB PDB, PubChem, Privacy, Terms — with the current page marked.
- **A breadcrumb that tracks real state**, not decoration: it names the loaded
  structure and the open section, and updates on both.
- **Section tabs are now a flat nav strip** rather than raised 3D tabs, which is
  what removes the tab-in-a-tab reading.
- **Content containers lost their two-tone bevels.** On the old tan window face
  a bevel read as a bevel; on white paper only its dark half is visible, so
  every box looked like it had a broken border. Containers now take a clean 1px
  rule. **Form controls keep their bevels** — a sunken input and a raised button
  are exactly what a form on such a site looked like.
- **Agency-style footer** — link row, a boxed research-use-only disclaimer, data
  source attribution, licence, and a "page generated" stamp — replacing the
  taskbar.
- **Added a skip-to-content link** as the first item in the tab order, and
  `aria-current` on the active section.

### Fixed
- Removed a dead `.raised`/`.sunken` helper pair that nothing referenced.


### Panels tab — the matrix engine gets a front end
- **Ligand × target heatmap** (`Panels` tab) — pick a curated system, run it, and
  read the matrix. Rows are ligands, columns are targets, and each cell is
  clickable for a full breakdown: raw score, ligand efficiency, z, rank in
  target, contact counts, top contact residues, which cofactors made it into
  the grid, and any measured ChEMBL activity for that pair.
- **The heatmap encodes z-scores, never raw scores.** This is the one thing the
  UI had to get right: a raw matrix read row-wise ranks pocket burial rather
  than preference, so the colour scale is bound to the server's per-target
  z-score and the caption says so.
- **Diverging scale, validated rather than eyeballed.** Polarity around z = 0
  (the panel's per-target average), two hues with a neutral midpoint. Each arm
  was checked as a single-hue ordinal ramp — monotone lightness, ≥ 0.06 step
  gaps, light end clearing the panel surface — and the poles verified
  CVD-separable (ΔE 24 protan / 31 normal vision).
- **Colour is never the only encoding.** Every cell carries its numeric z, the
  matrix is a real `<table>` with `scope`-d row/column headers and a caption,
  cells are keyboard-focusable and activate on Enter/Space, and measured
  activity is flagged with a glyph plus an `aria-label`, not a hue.
- **Failure and doubt are visible in the grid itself** — a target that could not
  be loaded gets a crimson cap and a struck-through header; a target whose
  declared cofactor is absent from its entry (so the site is scored as an empty
  cavity) gets an amber cap; failed cells are hatched; an unresolvable ligand
  is marked on its row label.
- **DRAFT curation banner** renders both when picking a system and on the result
  itself, since the result's verification block is the state it was actually
  computed under.
- **Live progress** — `submitJob` now surfaces the job's `progress` field, so a
  144-cell panel reports `Docking 47/144 (32.6%)` instead of sitting on
  "running" for minutes.

### CI
- **The front end is now syntax-checked on every push** (`node --check web/*.js`).
  `app.js` is the largest source file in the repo and was previously parsed only
  by the separate Playwright e2e workflow.


### Interface: "Luna Lab" theme
- **Full visual redesign** of `web/` as a Windows XP / late-90s scientific
  workstation: the app is now one window on a desktop — Luna-blue caption with
  working-looking caption glyphs, a menu bar, Explorer-bar task panels down the
  left rail, a real tab control, 3D-beveled buttons and sunken data wells,
  segmented "candy" progress bars, chunky period scrollbars, and a taskbar
  footer whose Start button carries the research-only warning.
- **Built entirely within the CSP.** `font-src 'self'` rules out webfonts, so
  the theme uses Tahoma / MS Sans Serif (period-correct anyway) with non-Windows
  fallbacks, and every bevel, gradient, caption glyph and grip is pure CSS or an
  inline `data:` SVG — no image assets and no new network origins.
- **Subpages inherit it for free** — privacy, terms and the API reference share
  `.topbar` + `.policy`, so they now render as floating XP dialogs on the
  desktop. They only needed a `subpage` body class to size the caption.
- **The 2D interaction diagram** was rethemed to the same palette (Luna-blue
  ligand node, pale-blue residue chips) and its label greys darkened for
  contrast against the lighter well.
- **Accessibility preserved and in two places improved** — visible focus rings
  kept (now amber, which reads on both the chrome and the panels),
  `prefers-reduced-motion` honoured, a `prefers-contrast: more` branch added,
  and two contrast defects fixed before they shipped: the taskbar Start button
  (3.2:1 → ~6:1) and the tray text (4.3:1 → ~4.8:1). The decorative menu bar is
  `aria-hidden`, so it is not announced as a menu that does nothing.
- **Narrow viewports drop the desktop metaphor** rather than shrink it: the
  window goes edge to edge, caption buttons hide, and the taskbar un-fixes so
  it cannot eat a small screen.
- Verified in headless Chromium across overview, interactions, docking, report,
  benchmark, a subpage and a 390px viewport, with **no uncaught JS errors** on
  the golden path. Status, tab, component and panel hooks are unchanged, so the
  existing e2e smoke test selectors still apply.


### Systems-level analysis: ligand x target interaction panels
- **Interaction panels** (`snaclex/panel.py`, job kind `panel`) — the workbench
  is no longer limited to one structure and one site. A panel docks **every
  ligand into every target** in a single run and returns a normalized matrix,
  so multi-molecule / multi-enzyme questions ("which of these metabolites
  engage which of these enzymes, and how selectively?") can be asked directly.
  One scoring grid is built per target and reused across all ligands, so the
  cost is M*N docks but only N grid builds.
- **Per-target normalization** — raw grid scores are *not* comparable across
  targets (a deeper pocket scores better for any ligand), so ligand efficiency
  is standardized within each target column (`z_target`, higher = better) and
  only those standardized values are compared across a ligand's row to give a
  `selectivity_gap`. The result states this explicitly; the previous batch
  screen never had to, because it only ever ranked within one site.
- **Measured activity alongside prediction** — each cell optionally carries the
  best matching ChEMBL activity for that (ligand, target) pair, so predictions
  sit next to experiment instead of standing alone.
- **Curated systems** (`snaclex/systems.py`, `snaclex/data/systems/*.json`,
  `GET /api/systems`) — versioned, reusable definitions of a piece of biology
  (targets + site hints + ligands + references) runnable as a panel from one
  id. Ships with a **catecholamine** system (12 targets x 12 ligands) covering
  synthesis, signalling and clearance, including two positive-control
  inhibitors so a panel's own rankings can be sanity-checked.
- **Verification is a first-class field, not a footnote** — a curated system
  asserts things a reader cannot check by running the code (that a PDB entry is
  the protein named, that the site ligand is really deposited, that the UniProt
  mapping holds). Every target therefore carries `verified: false` until
  `python -m snaclex.systems verify <id> --write` confirms it against RCSB, and
  **every API response and panel result carries the DRAFT flag** so an
  unverified curation cannot be mistaken for a checked one downstream. The
  shipped catecholamine system is a draft.
- **Job progress reporting** — `JobManager.submit(..., wants_progress=True)`
  injects a reporter, and `GET /api/jobs/{id}` now returns
  `{done, total, percent, label}`. A 144-cell panel runs for minutes; without
  this a client cannot tell a slow job from a wedged one.

### Fixed
- **Cofactors are no longer stripped from docking grids.**
  `docking.build_grid(..., extra_atoms=)` can now include heterocomponents in
  the rigid receptor. Previously only standard amino acids were scored, which
  silently removed FAD from a monoamine-oxidase site, SAM/Mg from
  catechol-O-methyltransferase and PLP from a decarboxylase — scoring
  cofactor-dependent sites as empty cavities. Curated targets declare their
  cofactors; any declared but absent from an entry is reported as a per-target
  warning rather than passing silently. The server's grid cache keys on the
  cofactor set, so a cofactor-free grid is never reused for a run that asked
  for one. Existing dock/screen/benchmark behaviour is unchanged (opt-in).

### Performance
- **Panels retain only the site shell of each receptor** (whole residues within
  `GRID_HALF + 12 A`), so a 16-target panel holds site neighbourhoods rather
  than 16 full assemblies. Proven exact by test: a docked pose profiles to
  identical contacts, counts and residues against the trimmed and untrimmed
  receptor. Residues are kept whole so aromatic-ring perception is unaffected.

### Documentation
- **`docs/systems-pharmacology.md`** — what a panel can and cannot answer, the
  three limits that bite hardest (cofactors, protonation, rigidity), and a
  five-rung capability ladder from here to dynamic pathway modelling, naming
  what each rung actually requires and where a structural workbench should stop
  and interoperate (SBML/COPASI) instead of reimplementing. Includes the
  validation strategy and where a spectral-prediction tool fits.

### Testing
- 204 offline tests (up from 176): the panel engine (shape, normalization,
  resilience, guards, progress, geometry), site-shell equivalence, cofactor
  handling, system validation/verification with injected RCSB fakes, and
  end-to-end panel jobs over HTTP.
- The API-contract test now **cross-checks every routed path against the
  documented contract**, so a new endpoint cannot ship undocumented (it
  previously only checked a hardcoded subset).


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
