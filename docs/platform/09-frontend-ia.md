# 09 — Frontend Information Architecture

Target: React + TypeScript + Mol*. Now: the existing vanilla-JS SPA in `web/`
grows a **Protein Report** shell alongside the preserved **Structural Analysis**
view. The IA is identical either way; only the rendering tech differs.

## 9.1 Top-level shell

```
┌ Global search bar ────────────────────────────────────────────────┐
│  [ protein / gene / accession / PDB / FASTA / chemical / variant ] │
│  detected type badge · organism filter · "paste a list" (batch)   │
└───────────────────────────────────────────────────────────────────┘
  Workspaces:  General  |  Immunology  |  Oncology  |  Genetics
  (workspaces are lenses over the same evidence model, not separate data)
```

On an ambiguous query the shell shows a **candidate picker** (organism,
accession, length, gene, review status) before opening a report (Stage-1 →
disambiguation, per [05](05-identifier-resolution.md)).

## 9.2 Protein Report — 13 tabs

Each tab reads from the `ProteinRecord`; every card exposes an **"i" provenance
affordance** linking to tab 13.

| # | Tab | Content | Key guardrails |
|---|---|---|---|
| 1 | **Overview** | preferred name, gene, organism, canonical accession, review status, length, structure availability, major functions, key domains, known ligands, **direct (A/B) vs predicted (D/E/F) chemical counts**, important variants, data freshness | Warning banners: ambiguous identity · low-confidence model · weak evidence |
| 2 | **Sequence & isoforms** | canonical sequence, isoform list + differences, checksum, sequence conflicts, length/MW | reviewed vs unreviewed badge; isoforms never collapsed |
| 3 | **Function & annotation** | function, catalytic activity, cofactors, localization, tissue specificity | evidence code + source per statement |
| 4 | **Domains & motifs** | domain/repeat/motif track (UniProt + InterPro), sequence feature ruler | source labelled (curated vs predicted) |
| 5 | **Structures** | hierarchy: experimental → homologous → predicted → none; per-structure coverage, %identity, chain map, missing residues, mutations, resolution/pLDDT, oligomeric state, ligands; **Mol* viewer** | AlphaFold badged "Predicted"; pLDDT coloring; docking-suitability flag |
| 6 | **Pockets & binding sites** | LIGSITE pockets (existing), volume/enclosure/druggability, conservation overlay | "geometric heuristic, not a confirmed site" |
| 7 | **Chemicals & ligands** | bound ligands (PDB), PubChem lookup, druglikeness, ligand→PubChem **match quality** (exact/parent/stereo/salt/...) | salts/parents/stereo not merged |
| 8 | **Protein interactions** | PPIs (IntAct/STRING/BioGRID, Phase 3), typed evidence | prediction vs curated separated |
| 9 | **Variants & mutations** | variant table, residue mapping, hotspots, domain disruption, WT-vs-mutant | clinical interpretation shows evidence + review status; never a diagnosis |
| 10 | **Pathways** | Reactome pathways (Phase 3) | source + evidence |
| 11 | **Disease & phenotype** | UniProt disease, ClinVar (Phase 3) | review status shown; association ≠ causation |
| 12 | **Literature** | PubMed references (Phase 3) | linked, not summarized as fact |
| 13 | **Evidence & provenance** | every evidence object; per-field source, version, retrieval date, method, params, limitations | the "why is this here?" tab |
| — | **Downloads** | JSON record, CSV tables, PDB/mmCIF, alignments, reproducible snapshot | checksums + tool/source versions |

## 9.3 Protein–chemical evidence panel (tab 7 / dedicated pair view)

Renders the A–F levels as **separate, badged cards** — never one score:

```
[A] Experimental structural   PDB 1IEP · chain A · ligand STI · contacts…       ● strongest
[B] Biochemical assay         AID 372 · Ki = 13 nM · active · Homo sapiens       ● measured
[C] Curated association       "imatinib inhibits ABL1" (source)                   ○ no value
[D] Homology-transferred      binds ABL1 homolog · 62% id · site 90% · limits…    ◐ inference
[E] Docking (this tool)       fit score −7.2 (lower=better) · seed 42 · params…   ◐ hypothesis, NOT affinity
[F] ML / similarity           target prediction · model vX                        ◔ weakest
Summary: direct A/B = 2 · predicted D/E/F = 3   (no boolean "interacts")
```

## 9.4 Workspaces (lenses)

- **Immunology** — cytokines/receptors, chemokines, checkpoints, Ig, antigen
  processing, **HLA/MHC allele-aware** (allele id never collapsed), epitopes
  (IEDB, gated), pathogen proteins, antibody–antigen structures, immune variants,
  immune pathways, immune-protein small molecules.
- **Oncology** — oncogenes/tumor suppressors, kinases, TFs, DNA-repair, apoptosis,
  cell-cycle; cancer variants + hotspots + residue mapping; pathway context;
  drug interactions (typed); resistance mutations; WT vs mutant structure
  coverage. Separates "protein associated with cancer" from "compound effective".
- **Genetics** — gene→transcript→protein, isoforms, variant consequences
  (missense/nonsense/frameshift/splice/start-loss/stop-loss), conservation,
  domain disruption, predicted structural effect, clinical interpretation *with
  evidence + review status + limitations*, ortholog comparison, phenotype links.

## 9.5 Search & comparison surfaces

Name/accession · sequence-similarity · motif · domain · organism filter · disease
& pathway filter · chemical substructure/similarity · by-bound-ligand · by-pocket
· by-variant-position · batch · family compare · ortholog compare · WT-vs-mutant ·
one-compound-many-targets · one-protein-many-compounds. Long searches show job
progress (async envelope).

## 9.6 Cross-cutting UI guardrails

- **Color/badge system** encodes evidence category everywhere: `experimental`
  (solid), `curated` (solid-muted), `calculated` (neutral), `homology` (hatched),
  `predicted`/`ml` (dashed). One legend, used site-wide.
- Predicted structures always say "Predicted (AlphaFold vN)"; never shown as
  experimental.
- Docking results always carry "fit score, not affinity; research-only".
- Every non-trivial number is click-through to its provenance.
- Empty/unavailable sources render as explicit "unavailable" states, not blanks.
