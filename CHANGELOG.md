# Changelog

All notable changes to SnaCleX are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

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
- **XSS defense-in-depth** — an `esc()` HTML-escaper for upstream-derived values
  interpolated into the DOM, complementing the CSP.

### Added
- **Privacy & Terms pages** (`web/privacy.html`, `web/terms.html`) linked from a
  new site footer: no accounts/cookies/trackers, what's collected, third-party
  data sources, retention, and the research-only disclaimer.
- **Test suite** (`tests/`) — 66 stdlib `unittest` cases covering the
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

_This covers Phase 0 (engineering hygiene) and Phase 1 (security & privacy
baseline) of the roadmap._

## [0.1.0]

- Initial public research tool: PDB structure loading (RCSB), atomic
  interaction profiling, LIGSITE pocket detection, Pfam-based conservation
  scoring, pure-Python Monte-Carlo docking, batch screening, PubChem/ChEMBL
  chemical lookup, and `.txt` session report export.
