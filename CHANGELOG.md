# Changelog

All notable changes to SnaCleX are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added — Motion module (Phase 4 of the new-modules brief)
- **Normal-mode analysis (ANM)** — a new **Motion** tab builds a Cα elastic
  network and solves the lowest-frequency normal modes (the large-scale
  hinge/breathing motions), reusing coordinates already loaded. `snaclex/motion.py`
  is pure analytic linear algebra — **not** molecular dynamics, no force field,
  no dependency. The 3N×3N Hessian's lowest internal modes are extracted by
  shift-invert subspace iteration with the 6 rigid-body modes projected out
  analytically (a stdlib Cholesky + Jacobi eigensolver, cross-checked against a
  brute-force full solve). Large structures are coarse-grained to ≤60 beads to
  stay synchronous (documented cap). `/api/motion`, cached.
- **Viewer animation** — mode selector + amplitude slider; "Animate" oscillates
  the structure along the selected mode in the existing 3Dmol viewer via
  multi-frame models.
- **Contact-strain tie-in (the stronger feature)** — under the selected mode and
  amplitude, the module re-evaluates the **Interactions** contacts of a bound
  molecule and ranks which strain/break most, turning "here is a wiggle" into
  "here is the interaction weak point of this motion."
- **MD explicitly deferred**; report section, JSON export, and
  `provenance.motion_methods()` (method, cap, relative-frequency caveat) added.

### Added — HLA / MHC module (Phase 3 of the new-modules brief)
- **MHC-fold detection on load** — structures are classified as MHC **class I**
  (polymorphic heavy chain + β2-microglobulin) or **class II** (α + β chains) by
  chain alignment to reference sequences, surfacing a new **HLA / MHC** tab and an
  "MHC class N detected" badge. `snaclex/hla.py`, dependency-free.
- **Reference-annotated groove** (NOT blind LIGSITE) — the peptide-binding groove
  is annotated from a fixed reference position set (class I pockets **A–F**;
  Saper/Bjorkman/Wiley 1991, Madden 1995), mapped onto the loaded allele's
  numbering by sequence alignment and highlighted with a new `HLA groove pockets`
  color mode + legend (reusing the residue-coloring path). Class II is
  best-effort (β-chain).
- **Allele-specific groove comparison** — each groove position is compared to a
  class-I reference (A*02:01); allele-specific residues are reported per pocket
  with charge / hydrophobicity / size deltas (reusing the Pockets residue-property
  sets) and drawn as sticks in 3D.
- **HLA↔drug hypersensitivity lookup** — a curated table (B*57:01↔abacavir,
  B*15:02 & A*31:01↔carbamazepine, B*58:01↔allopurinol, plus dapsone/vancomycin),
  each row cited to primary literature. Hits surface a **"Look up drug →"** button
  that pulls the drug up in the existing Chemical/PubChem viewer. A lookup, not a
  prediction; **RESEARCH ONLY** banner prominent (clinically adjacent). AFND /
  IPD-IMGT/HLA credited as the population-frequency / sequence sources.
- **Peptide-into-groove docking deferred** behind an experimental flag (kept out
  of the synchronous compute budget), clearly labeled.
- Pockets tab now notes the groove is reference-annotated (don't blind-detect it);
  report section, JSON export, and `provenance.hla_methods()` added.

### Added — Genome variant bridge (Phase 2 of the new-modules brief)
- **ClinVar / gnomAD overlay** — a new **Variants** tab resolves the loaded
  structure to its UniProt accession, pulls missense variants from the EMBL-EBI
  Proteins API (aggregating **ClinVar** clinical significance and **gnomAD**
  population allele frequencies), and maps each variant's UniProt position onto
  the structure's residue numbering by aligning each chain to the UniProt
  sequence (reusing the Evolution NW aligner). `snaclex/variants.py` +
  `/api/variants`, cached per structure. This is a genotype→structure *bridge* —
  no genome browser, ideogram, coordinate system, or reads.
- **Two coloring modes** reusing the residue-coloring path — a new
  `Variants (ClinVar/gnomAD)` viewer mode with a **pathogenicity ⇄ frequency**
  toggle; clicking a colored residue shows its substitutions, ClinVar
  significance and gnomAD AF.
- **Interpretation stack** — each variant residue is cross-annotated (client
  side) with whether it also falls in a detected **pocket**, on a **conserved**
  residue, or at an **interaction contact** already computed this session.
- **Gene/locus card** — gene symbol, chromosome and cytogenetic band via NCBI
  E-utilities (best-effort, info-only, degrades gracefully).
- **RESEARCH-ONLY** banner is prominent on the tab (clinically adjacent data),
  with report section, JSON export, and `provenance.variant_methods()`
  (source + retrieval date + limitations).

### Added — Antibody layer (Phase 1 of the new-modules brief)
- **Antibody detection on load** — every loaded/uploaded structure is scanned
  for immunoglobulin variable domains. Each protein chain is typed (VH, VL,
  heavy/light constant, or "other" candidate antigen) by aligning its N-terminal
  window to germline framework consensus references and requiring the two
  invariant Ig-domain cysteines. A small `snaclex/antibody.py` module does this
  dependency-free, reusing the existing Needleman–Wunsch aligner; it runs inline
  in `/api/analyze` and `/api/upload` (single-digit ms) and drives a new
  **Antibody** tab plus an "antibody detected" badge.
- **CDR highlighting** — CDR1/CDR2/CDR3 are delimited (IMGT) and recolored in the
  3D viewer through the **same residue-coloring path** the Evolution module uses
  (new `CDR loops (antibody)` color mode + legend).
- **Developability liability scan** — sequence-motif flags (N-glycosylation
  `N-X-S/T`, deamidation `N-G/N-S`, isomerization `D-G`, unpaired cysteine), each
  marked in-CDR (higher concern) vs framework.
- **Paratope–epitope interface** — when an antibody chain and a candidate antigen
  coexist in one file, `/api/antibody_interface` computes the protein–protein
  contact set (paratope vs epitope residues) and a coarse Shrake–Rupley buried
  surface area, as a conditional branch of the interaction logic — instead of
  routing a flat paratope through LIGSITE.
- **Interpretation caveats** — the Evolution tab now warns that Pfam alignments
  underrate hypervariable CDRs for antibodies, and the Pockets tab marks its
  output low-confidence when an antibody is loaded (LIGSITE mislabels paratopes).
- **Report + provenance** — an antibody section is added to the compiled report
  and the JSON export, with `provenance.antibody_methods()` (numbering scheme,
  method, limitations, research-only disclaimer).

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
