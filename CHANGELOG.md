# Changelog

All notable changes to SnaCleX are documented here. This project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **Test suite** (`tests/`) — 55 stdlib `unittest` cases covering the
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

_This work is Phase 0 (engineering hygiene) of the roadmap._

## [0.1.0]

- Initial public research tool: PDB structure loading (RCSB), atomic
  interaction profiling, LIGSITE pocket detection, Pfam-based conservation
  scoring, pure-Python Monte-Carlo docking, batch screening, PubChem/ChEMBL
  chemical lookup, and `.txt` session report export.
