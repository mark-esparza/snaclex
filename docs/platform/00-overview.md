# SnaCleX Platform Expansion — Design Set

This directory contains the design deliverables for expanding SnaCleX from a
PDB-centred docking workbench into a **sequence-first, structure-aware protein &
chemical interaction platform**. It is the paper trail requested by the task:
produce the design before the code, keep the existing structural tooling intact,
and never overstate evidence.

## Read in this order

| # | Document | Deliverable |
|---|----------|-------------|
| 01 | [Product requirements](01-prd.md) | PRD |
| 02 | [Architecture & data flow](02-architecture.md) | System architecture + data-flow diagrams |
| 03 | [Source integration table](03-source-integration.md) | Source-by-source integration table |
| 04 | [Normalized data model](04-data-model.md) | Protein Record, Evidence object, knowledge graph |
| 05 | [Identifier resolution](05-identifier-resolution.md) | Identifier-resolution strategy |
| 06 | [Evidence ranking](06-evidence-ranking.md) | Evidence-ranking specification (levels A–F) |
| 07 | [API definitions](07-api.md) | Versioned API endpoints |
| 08 | [Backend structure](08-backend-structure.md) | Backend repository structure |
| 09 | [Frontend IA](09-frontend-ia.md) | Frontend information architecture (13 report tabs) |
| 10 | [Phased plan](10-phased-plan.md) | Phased implementation plan + vertical slice |
| 11 | [Researcher workflows](11-workflows.md) | Example researcher workflows |
| 12 | [Deployment](12-deployment.md) | Deployment instructions |
| 13 | [Validation plan](13-validation-plan.md) | Scientific validation plan + example cases |

The **working vertical slice** that accompanies these documents lives in the
`snaclex/` package (new modules) and `server.py` (new endpoints); see
[10-phased-plan.md](10-phased-plan.md) §"Vertical slice — what shipped".

## Core principles (non-negotiable)

1. **Sequence first.** A `ProteinRecord` exists independently of whether any 3D
   structure is available. Every analysis path degrades to sequence-only rather
   than failing.
2. **Evidence is typed and never flattened.** Direct experimental structural
   evidence, biochemical assays, curated associations, homology transfer,
   docking, and ML predictions are six *distinct* levels (A–F). A docking score
   is never presented as a measured affinity; an association is never presented
   as proof of binding.
3. **Provenance on every field.** Source, source identifier, retrieval date,
   source version, evidence category, and limitations travel with each value.
4. **Identity is resolved, not guessed.** Proteins are never merged on gene
   symbol or name alone — accession, taxonomy, sequence checksum, isoform, and
   cross-reference mappings decide identity. One-to-many and many-to-many links
   are preserved.
5. **Preserve the existing capability.** The current RCSB structure loading,
   interaction profiling, pocket detection, conservation, and docking become a
   modular **Structural Analysis** component; they are refactored around, not
   removed.
6. **Zero-dependency default, opt-in heavyweight track.** See "Stack decision".

## Stack decision (explicit assumption)

The task's recommended default stack is FastAPI / Pydantic / PostgreSQL /
SQLAlchemy / Redis / Celery / S3 / React+TypeScript / Mol*. The task also says:
*"Use different technologies only when there is a documented technical
advantage."*

SnaCleX has a **documented, deliberate zero-dependency design**
(`requirements.txt` empty by design, `render.yaml` has no build step, README:
"Pure Python standard library — no pip installs required"; see `ROADMAP.md`
§"Guiding constraint"). Discarding that to bolt on Postgres/Redis/Celery/React
would remove the project's main operational virtue (deploys anywhere, no supply
chain, fully reproducible offline tests).

**Decision:** the *target production architecture* is specified against the
recommended stack (see [02](02-architecture.md), [08](08-backend-structure.md)),
because it is the right destination at scale. The **vertical slice ships in the
existing stdlib idiom** so it runs today with no new dependencies and its tests
stay fully offline. Each service is defined behind an interface so the stdlib
implementation can be swapped for the FastAPI/Postgres/Celery one without
touching callers. This mirrors the repo's existing "stdlib-first / optional
heavyweight track" split.

## Assumptions (called out explicitly, per the task)

- **A1.** Public API shapes below reflect the current documented endpoints of
  UniProtKB, UniParc, RCSB, NCBI E-utilities, AlphaFold DB, and PubChem as of
  2026-07. Any field marked **`⚠ VERIFY`** in the source table or code comments
  has *not* been round-tripped against a live response in this environment
  (offline CI) and must be confirmed before relying on it in production.
- **A2.** No redistribution of source data is assumed. Adapters fetch on demand
  and cache locally; bulk ingestion/mirroring is a later, license-gated step
  (see [03](03-source-integration.md) "License / attribution").
- **A3.** On-demand structure *prediction* (ESMFold/ColabFold/local AlphaFold)
  is designed as an optional, env-gated interface only. It is never required and
  is not part of the vertical slice.
- **A4.** No clinical, diagnostic, or therapeutic claims are produced anywhere.
  Variant and disease views display originating evidence, review status, and
  limitations; they never interpret a variant clinically.
- **A5.** "Homolog-transferred" and "ML-predicted" evidence are always labelled
  as inference, with the transfer basis (sequence identity, aligned site
  identity) shown.

## What is intentionally *not* built yet

Only Phase 1 (unified sequence-first record + evidence separation) has a working
slice. Phases 2–4 (AlphaFold-confidence-aware docking, InterPro, workspaces,
knowledge-graph exploration, on-demand prediction) are specified but scaffolded,
not implemented, to honour the task's "do not generate a large amount of
disconnected code" instruction.
