# 13 — Scientific Validation Plan

Two goals: (1) confirm each source integration returns correct, current data, and
(2) confirm the platform **never overstates evidence** across the 12 required
case categories. Validation is layered: unit → integration → scientific-rubric →
regression.

## 13.1 Validation layers

| Layer | Scope | How (matches existing conventions) |
|---|---|---|
| **Unit (offline)** | pure-compute + parse helpers | stdlib `unittest`, recorded fixtures; no network. Existing `tests/` pattern. |
| **Adapter contract (offline)** | each adapter's parse over a recorded real response | fixture JSON captured once from the live API, checked into `tests/fixtures/`; assert normalized fields + provenance. |
| **Integration (in-process)** | new endpoints end-to-end | boot `ThreadingHTTPServer`, mock the HTTP seam, assert response shape + guardrails (as in `test_server_integration.py`). |
| **Live smoke (opt-in, CI-gated)** | real API reachability + field presence | a `@network` test suite run manually / nightly; flags `⚠ VERIFY` fields that drifted. |
| **Scientific rubric (manual)** | the 12 case categories below | reviewer checklist: useful output? evidence not overstated? provenance present? |
| **Regression** | golden outputs | snapshot the report JSON for canonical cases; diff on change. |

## 13.2 Scientific safeguards checklist (run per release)

For every case, a reviewer confirms:

- [ ] Experimental facts are visually/structurally distinct from computational
      predictions.
- [ ] Source-specific evidence is shown for each major statement.
- [ ] Isoform and organism context preserved (nothing collapsed on name/gene).
- [ ] Structure coverage and sequence identity reported.
- [ ] Model-confidence metrics (pLDDT/PAE) exposed for predicted models.
- [ ] No causal claim derived from association data.
- [ ] No docking score rendered as biological activity / affinity.
- [ ] No clinical recommendation anywhere.
- [ ] No in-silico result described as "validation".
- [ ] Conflicting records reported, not silently resolved.
- [ ] Negative / inactive assay results displayed when relevant.
- [ ] Assay context + units preserved.
- [ ] Database + software versions tracked; analysis reproducible from snapshot.

## 13.3 Required case categories (from the task) + acceptance

| # | Category | Example probe | Accept when… |
|---|---|---|---|
| 1 | Many experimental structures + ligands | HIV-1 protease (`1HSG`), carbonic anhydrase II (`1CA2`) | structures + bound ligands listed; existing interaction/pocket output intact |
| 2 | Predicted structure, no experimental | a Swiss-Prot protein with only `AF-…-F1` | AlphaFold model shown, badged Predicted, pLDDT bands present; no PDB claimed |
| 3 | No usable structure | disordered/uncharacterized TrEMBL entry | sequence-only analysis returns; docking gated off; banners raised |
| 4 | Multiple isoforms | BRCA1 (P38398) | isoforms listed separately, not merged |
| 5 | Clinically studied variant | `BRAF V600E` (P15056) | residue mapped (seq + structure coords); interpretation shows evidence + review status, no diagnosis |
| 6 | Human immune protein | PD-1 (`PDCD1`), an HLA allele | immunology lens; HLA allele identity preserved |
| 7 | Cancer-associated kinase | ABL1 (`P00519`), `T315I` | kinase domain + hotspot; WT vs mutant coverage; construct kept in evidence |
| 8 | Microbial enzyme | bacterial DHFR via `WP_` RefSeq | organism/strain preserved; RefSeq→UniParc→UniProt bridge resolves |
| 9 | Poorly characterized protein | a hypothetical/DUF protein | useful sequence-level output; low-annotation state honest |
| 10 | Chemical with direct assay evidence | imatinib (CID 5291) × ABL1 | Level B assay values incl. inactives; not flattened |
| 11 | Chemical with only homolog evidence | tool compound × uncharacterized homolog | Level D only, with transfer basis; no A/B claimed |
| 12 | Protein–chemical pair, no prior evidence | arbitrary novel pair | A/B/C/D empty reported honestly; Level E offered, labelled hypothesis |

## 13.4 Known canonical anchors (reuse existing golden tests)

The repo already validates the structural core against real cases; keep them as
regression anchors and extend:

- Pocket detection recovers the true ligand site as top pocket for `1HSG`/`1CA2`
  (existing).
- Redocking benzamidine into trypsin `3PTB` reproduces the S1 pose (Asp189 salt
  bridge, ~2 Å redock RMSD) (existing).
- **New:** raw-FASTA-of-P38398 resolves (checksum bridge) to UniProt P38398 (
  adapter contract test with recorded UniParc `checksum:` response).
- **New:** an accession with no PDB but an AlphaFold model yields a
  `structures.predicted[0]` and empty `structures.experimental` (mocked).
- **New:** an evidence assembly test asserting a docking item is Level E with
  `units == "score"` and can never carry `Ki/nM` (invariant 6.4-1).

## 13.5 `⚠ VERIFY` retirement process

Every field marked `⚠ VERIFY` in [03](03-source-integration.md) / code comments
must be confirmed against a live response before it is depended on:

1. Capture a live response into `tests/fixtures/<source>_<case>.json`.
2. Add an adapter-contract test asserting the field path + type.
3. Remove the `⚠ VERIFY` marker in the same commit.

Until retired, code treats such fields as best-effort: absence → `null` +
provenance note, never a fabricated value.
