# New-modules implementation notes

How the four planned branches (Antibody, Genome-variant bridge, HLA/MHC, Motion)
map onto SnaCleX's actual architecture. Written before coding, per the brief's
"FIRST STEPS". Where the brief and the codebase disagreed, the codebase wins —
those points are flagged.

## How the current app is wired (the spine every new module reuses)

- **Backend**: Python-stdlib-only. `server.py` is a `ThreadingHTTPServer` with a
  hand-rolled router (`_route`). Each analysis lives in a `snaclex/<name>.py`
  module exposing a pure function over a parsed `Structure`. The server adds a
  cached `_get_<name>(pdb_id)` helper and a `_api_<name>(qs)` handler.
- **Structure model**: `snaclex/pdbparse.py` → `Structure(atoms, protein_atoms,
  components, chains)`; `Atom` / `Component` dataclasses. PDB **and** mmCIF are
  parsed to the same model.
- **Sequence / alignment utilities**: `evolution.AA3TO1`, `structure_sequence`,
  and the Needleman–Wunsch aligner `evolution._nw_align` are reused directly.
- **Contact geometry**: `interactions._Grid` + `_dist` + the cutoffs (`HB_MAX`,
  …). Protein–protein interface analysis is a conditional branch here, not a new
  engine.
- **Residue coloring in 3D**: `web/app.js` `rebuildScene()` styles residues with
  `viewer.setStyle({chain, resi:[…]}, {cartoon:{color}})`, grouped by shade
  (`applyConservationColors`). Every new "color by X" reuses this path and adds
  an `<option>` to `#colorMode`.
- **Provenance**: `snaclex/provenance.py` returns static method blocks; the
  report (`app.js` `buildReportSections`) and the JSON/TXT exports embed them.
- **Report**: `buildReportSections()` aggregates `state.*`; `compileReport()`
  renders; exports carry provenance. New modules push a section here.

## Phase 1 — Antibody layer (this change)

- **`snaclex/antibody.py`** (new): lightweight, dependency-free.
  - Chain typing by aligning each chain's N-terminal sequence to germline
    framework consensus references (VH, Vκ, Vλ, and constant CH1/CL) with the
    existing `_nw_align`. Typing keys off framework identity **plus** the two
    invariant Ig-domain cysteines — specific enough not to fire on every
    Ig-superfamily fold. Chains become `VH` / `VL` / `heavy_constant` /
    `light_constant` / `other` (candidate antigen).
  - CDR delimitation (IMGT by default) is read off the reference columns, with
    CDR3's variable length captured by the "last framework segment seen" walk of
    the alignment and snapped to the J-motif (W/F-G-x-G).
  - Liability scan: pure regex over the VH/VL sequence (N-glycosylation
    `N[^P][ST]`, deamidation `N[GS]`, isomerization `DG`, unpaired Cys), each
    flagged in-CDR vs framework.
- **`interactions.profile_interface`** (new fn, same module): reuses `_Grid` to
  collect inter-chain heavy-atom contacts → paratope vs epitope residue sets,
  with a coarse Shrake–Rupley buried-surface-area estimate over interface
  residues only (bounded compute). This is the "antibody + antigen in one file"
  branch that replaces routing a flat paratope through LIGSITE.
- Detection runs inline in `/api/analyze` (cheap: NW on a few short chains), so
  the badge appears on load. Interface analysis is an on-demand endpoint
  (`/api/antibody_interface`) like pockets/evolution.
- **Frontend**: new "Antibody" tab, per-chain badges, a `cdr` color mode +
  legend, liability table, the Evolution/Pfam CDR caveat (1.5), and a
  low-confidence note on Pockets when an antibody is loaded (1.3 design note).
- **Provenance**: `provenance.antibody_methods()`; report + exports updated.

Deferred (backlog, per brief): germline/humanness, Chothia canonical classes,
epitope binning.

## Phase 2 — Genome variant bridge (shipped)

- **`snaclex/variants.py`** (new): resolves PDB→UniProt via the existing
  `rcsb.fetch_uniprot_accessions`, then pulls missense variants from the EMBL-EBI
  Proteins API variation endpoint — which aggregates **ClinVar** clinical
  significance and **gnomAD** allele frequencies *and* ships the UniProt sequence.
  Each chain is NW-aligned to that sequence (reusing `evolution._nw_align`) so
  UniProt positions land on structure residues (works across author numbering and
  homomers). Cached per structure in `server._get_variants`; `/api/variants`.
- **Frontend**: a **Variants** tab, a `variant` color mode with a
  pathogenicity⇄frequency toggle (both reuse the `setStyle` residue path),
  variant details on residue click, a gene/locus card (NCBI E-utilities,
  best-effort), and client-side cross-annotation against pockets / evolution /
  interactions. RESEARCH-ONLY banner is prominent (clinically adjacent).
- Same EBI host the Evolution module already depends on, so no new
  connectivity assumptions; all fetches degrade gracefully to "data unavailable".

## Phase 3 — HLA / MHC module (shipped)

- **`snaclex/hla.py`** (new): MHC-fold detection (class I heavy + β2m; class II
  α/β) by chain alignment to reference sequences — runs inline on load like the
  antibody layer. Groove annotated from a fixed reference pocket-lining set
  (class I A–F), mapped onto the structure by NW alignment to an A*02:01
  reference (so it is reference-annotated, not blind LIGSITE). Allele-specific
  groove comparison reuses the Pockets residue-property sets. Curated HLA↔drug
  hypersensitivity table with primary-literature citations. `/api/hla` (cheap,
  synchronous). Peptide docking deferred behind `PEPTIDE_DOCKING_AVAILABLE`.
- **Frontend**: HLA/MHC tab (allele input), `groove` color mode + pocket legend,
  groove/pocket + allele-difference tables, drug hits with a "Look up drug"
  button into the existing Chemical viewer, MHC badge, Pockets note, report +
  export. RESEARCH-ONLY banner prominent (clinically adjacent).
- IPD-IMGT/HLA + AFND credited as sequence/frequency sources. EBI/live allele
  sequence fetch is not required — the loaded structure supplies the allele's
  actual groove residues.

## Phase 4 (planned, not yet built)

- **Motion (NMA/ANM)**: analytic ANM on Cα coordinates (numpy-free linear
  algebra), mode selector + amplitude slider animating the viewer, and a
  contact-strain readout tied back to Interactions.
