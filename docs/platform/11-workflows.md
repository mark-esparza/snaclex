# 11 — Example Researcher Workflows

Concrete end-to-end walkthroughs against the platform. Each maps to the pipeline
stages and shows the scientific guardrails in action.

## W1 — "I have a gene symbol, no structure in mind" (BRCA1)

1. Enter `BRCA1`. Stage-1 classifies `gene_symbol`; candidate picker shows
   *Homo sapiens P38398 (reviewed, 1863 aa)* first, other organisms listed.
2. Pick human. Stage-2 builds the xref graph
   (UniProt→UniParc→RefSeq `NP_009225`→PDB entities→`AF-P38398-F1`).
3. Report opens: **Overview** shows reviewed status, functions, BRCT domains,
   *partial* experimental coverage (only some domains crystallized) → **Structures**
   tab lists experimental fragments + full-length AlphaFold model (badged
   Predicted), plus a "no single full-length experimental structure" note.
4. **Sequence & isoforms** lists isoforms without collapsing them.
5. Researcher exports the JSON snapshot (provenance-stamped).

*Demonstrates:* sequence-first record, partial structure coverage, predicted
fallback, isoform preservation.

## W2 — "Protein with no usable structure" (a disordered/uncharacterized protein)

1. Enter a TrEMBL accession. Overview raises **"unreviewed record"** banner.
2. **Structures** tab: no PDB; AlphaFold present but large low-confidence
   regions → `docking_suitable: false`, "sequence-only; docking not run".
3. **Sequence analysis** still returns composition, MW, pI, charge@pH, GRAVY,
   low-complexity — useful output with no structure.

*Demonstrates:* graceful degradation, docking gate on low confidence, "useful
without overstating".

## W3 — "Protein–chemical pair with strong evidence" (ABL1 + imatinib)

1. Enter pair `ABL1 + imatinib`. Compound normalizes to PubChem CID (parent);
   `chemical_form.matched_as` recorded.
2. **Evidence panel** renders separate cards:
   - **A** — imatinib bound to ABL1 in PDB (e.g. `1IEP`), chain/ligand/contacts.
   - **B** — PubChem BioAssay Ki/IC50 values, including any inactive outcomes.
   - **E** — optional docking into the experimental site → *fit score*, seed,
     params; explicitly "not an affinity".
3. Overview counts **"direct A/B = 2, predicted D/E = 1"** — no boolean
   "interacts".

*Demonstrates:* A–F separation, docking≠affinity, negatives retained.

## W4 — "Compound with only homolog evidence" (tool compound vs an uncharacterized kinase)

1. Enter pair; no PDB or assay for *this* protein.
2. Evidence gathers **Level D**: compound binds a 62%-identity homolog; card
   shows source protein, sequence identity, aligned binding-site identity, and
   transfer limitations. No A/B claimed.

*Demonstrates:* honest homology transfer with its basis shown.

## W5 — "Protein–chemical pair with no prior evidence"

1. Enter pair. A/B/C/D all empty → panel shows explicit "no direct/curated/
   homolog evidence found".
2. If the protein has a suitable structure, **Level E docking** is *offered* (not
   auto-run) as a hypothesis-generation step, clearly labelled.

*Demonstrates:* absence reported honestly; docking as hypothesis only.

## W6 — "A clinically studied variant" (BRAF V600E)

1. Enter `BRAF V600E`. Stage-1 → `variant`; parsed to `BRAF`(P15056) residue 600
   V→E.
2. **Variants** tab: residue mapped in sequence and (if structure covers it) in
   structure coordinates; domain context (kinase domain); WT-vs-mutant structure
   coverage; any curated/clinical interpretation shown **with review status and
   limitations** — never a clinical recommendation.
3. **Oncology workspace** lens adds hotspot context and typed drug associations.

*Demonstrates:* variant→residue mapping, dual coordinates, no clinical call.

## W7 — "Cancer-associated kinase, WT vs mutant" (ABL1 T315I)

1. Enter `ABL1 T315I`. Compare WT vs mutant structure coverage.
2. Evidence for imatinib shows **Level B** with `protein_construct: "T315I"` kept
   distinct from wild-type assays; resistance context surfaced in the oncology
   lens.

*Demonstrates:* construct/isoform context preserved in evidence.

## W8 — "Human immune protein" (PD-1 / HLA)

1. Enter `PDCD1` → immunology workspace; checkpoint context.
2. For an HLA query (`HLA-A*02:01`), the record identity is the **allele**, not
   the gene — alleles are never collapsed.

*Demonstrates:* allele-aware immunology handling.

## W9 — "Microbial enzyme" (e.g. a bacterial dihydrofolate reductase)

1. Enter organism+name or a WP_ RefSeq accession. Organism/strain preserved;
   RefSeq→checksum→UniParc→UniProt bridge resolves identity.

*Demonstrates:* non-human, RefSeq-first resolution.

## W10 — "Batch" (a list of proteins)

1. Paste a newline list. Each resolves independently (async job); results table
   links to individual reports. Family/ortholog comparison available.

*Demonstrates:* batch + comparison surfaces.

Each workflow above is realizable with the shipped slice at the API level for
W1–W5/W9 (protein + evidence endpoints); W6–W8/W10 use Phase-2/3 surfaces that
are specified and scaffolded.
