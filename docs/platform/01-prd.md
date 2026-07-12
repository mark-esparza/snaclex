# 01 — Product Requirements Document

## 1. Problem

SnaCleX today answers *"what does this chemical do inside this crystal
structure?"* It requires an experimentally determined PDB entry as the entry
point. Researchers routinely start from a **protein sequence or identifier**
(a gene symbol, a UniProt accession, a raw FASTA, a mutation like `BRAF V600E`)
for proteins that may have **no experimental structure at all**. For those, the
current tool has no answer.

## 2. Vision

A **sequence-first, structure-aware** platform that accepts any publicly
cataloged protein or chemical identifier and returns a provenance-rich report —
sequence analysis always, structural analysis when a structure (experimental,
homologous, or predicted) exists, and **typed** protein–chemical interaction
evidence that never conflates measurement with prediction.

## 3. Goals

- **G1.** Accept any of: protein name, gene symbol, UniProt accession, UniParc
  UPI, NCBI Protein / RefSeq accession, PDB ID, organism+name, raw FASTA, a list
  of proteins, a PubChem CID / name / SMILES / InChI / InChIKey, a
  protein–chemical pair, or a variant (`BRAF V600E`, accession+substitution).
- **G2.** Build a unified `ProteinRecord` that exists with or without structure.
- **G3.** Resolve identity across UniProtKB / UniParc / NCBI / RCSB / AlphaFold /
  PubChem without collapsing isoforms, paralogs, orthologs, chains, or fragments.
- **G4.** Produce a structure-availability hierarchy (experimental → homologous →
  predicted → sequence-only) with coverage, chain mapping, sequence identity, and
  confidence.
- **G5.** Produce a **6-level evidence hierarchy (A–F)** for every
  protein–chemical pair and never merge levels.
- **G6.** Preserve and modularize the existing structural tooling.
- **G7.** Make every major statement inspectable: source, record, method, params,
  retrieval date, limitations.

## 4. Non-goals

- Not a clinical decision tool. No diagnosis, prognosis, or therapy advice.
- Not an affinity predictor. Docking scores are fit heuristics, not kcal/mol.
- Not a data mirror. No wholesale redistribution of source databases.
- On-demand structure prediction is optional and out of the initial scope.

## 5. Personas

| Persona | Needs |
|---|---|
| **Structural biologist** | Existing PDB workflow, pockets, interactions, docking — unchanged. |
| **Molecular / cell biologist** | Function, domains, localization, interactions for a gene with maybe no structure. |
| **Cancer / immunology researcher** | Variant→residue mapping, hotspots, allele-aware HLA, checkpoint proteins, drug associations with *typed* evidence. |
| **Cheminformatician / pharmacologist** | Compound normalization (parent/salt/stereo), bioassay evidence, target comparison. |
| **Geneticist** | Gene→transcript→protein, variant consequences, conservation, clinical interpretations *with evidence shown*. |

## 6. Functional requirements (by pipeline stage)

- **FR-1 Query interpretation.** Classify input type; present candidate matches
  (organism, accession, length, gene, review status) on ambiguity.
- **FR-2 Identifier normalization.** Cross-reference graph across UniProt / UPI /
  NCBI / RefSeq / Gene / PDB entity+chain / AlphaFold / PubChem, preserving
  one-to-many and many-to-many edges.
- **FR-3 Sequence analysis (always).** Length, MW, composition, pI, charge at a
  chosen pH, hydrophobicity/GRAVY, low-complexity, plus retrieved features
  (signal peptide, TM, domains, active/binding sites, PTMs, variants). Every
  result labelled calculated / curated / homology-inferred / ML-predicted.
- **FR-4 Structure availability.** Experimental (PDB) → homologous → predicted
  (AlphaFold) → on-demand (optional) → sequence-only, with per-structure
  coverage, chain map, %identity, mutations, missing residues, confidence,
  oligomeric state, bound ligands, docking suitability.
- **FR-5 Structural analysis.** The existing module, plus confidence-aware
  handling of predicted models; residue numbering retained in both source and
  canonical coordinates.
- **FR-6 Protein–chemical evidence.** Gather A–F evidence before docking;
  present each level distinctly.
- **FR-7 Chemical normalization.** PubChem CID / InChIKey / parent / stereo /
  protonation / salt / tautomer; classify PDB-ligand→PubChem match quality.
- **FR-8 Report.** 13 tabbed sections (see [09](09-frontend-ia.md)) with a
  provenance tab and downloadable, reproducible snapshot.

## 7. Non-functional requirements

- **NFR-1 Reproducibility.** Every analysis records tool + source versions and
  retrieval timestamps; a snapshot can be re-opened without re-querying.
- **NFR-2 Resilience.** Any single source failing (or being license-restricted)
  degrades gracefully; the core never hard-depends on an optional adapter.
- **NFR-3 Rate/limits.** Per-source rate limiting, retry with backoff,
  conditional requests where supported, optional API keys (NCBI, others).
- **NFR-4 Separation of paths.** Ingestion / cache-refresh / interactive query
  are separate; a report open must not fan out live calls to every source.
- **NFR-5 Scientific safety.** Experimental vs computational always
  distinguished; no causal claims from association; conflicting records reported,
  not silently resolved; negative/inactive assays shown when relevant.
- **NFR-6 Zero-dependency default.** Core runs on the Python standard library;
  heavyweight integrations are opt-in and env-gated.

## 8. Success metrics

- A researcher can get a useful report for a protein with **no** experimental
  structure (AlphaFold + sequence-only) — the headline capability gap closed.
- For the 12 validation cases ([13](13-validation-plan.md)), the platform
  returns useful output **without overstating evidence** (manual rubric).
- Zero instances of a docking score rendered as an affinity, or an association
  rendered as proven binding (enforced by copy review + unit tests on labels).

## 9. Constraints & risks

- **C1.** Public APIs change field names (already seen: PubChem SMILES field
  rename handled in `pubchem.py`). Adapters must degrade gracefully and be
  version-pinned where possible.
- **C2.** Homology transfer and ML prediction are easy to over-trust; the UI must
  visually separate them (color/badge) and always show the transfer basis.
- **C3.** Licensing differs per source (see [03](03-source-integration.md)); some
  enrichment sources cannot be redistributed and must be fetched live/attributed.
- **C4.** Isoform/allele collapsing is a correctness risk (esp. HLA/MHC); identity
  resolution ([05](05-identifier-resolution.md)) is the mitigation.
