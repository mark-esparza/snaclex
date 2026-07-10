# 05 — Identifier-Resolution Strategy

Covers pipeline **Stage 1 (query interpretation)** and **Stage 2 (identifier
normalization)**. Implemented in the slice as `snaclex/idresolve.py`.

## 5.1 Query interpretation (classification)

Input is classified by cheap, ordered pattern rules *before* any network call.
Ambiguous inputs return **candidates**, never a silent pick.

| Rule (checked in order) | Pattern | Classified as |
|---|---|---|
| FASTA | starts with `>` or is a run of ≥20 valid AA letters | `raw_sequence` |
| PDB id | `^[0-9][A-Za-z0-9]{3}$` (and not a UniProt shape) | `pdb_id` |
| UniProt acc | `^[OPQ][0-9][A-Z0-9]{3}[0-9]$` or `^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$` (opt. `-\d+` isoform) | `uniprot_accession` |
| UniParc UPI | `^UPI[0-9A-F]{10}$` | `uniparc_upi` |
| RefSeq protein | `^(NP|XP|WP|YP|AP)_\d+(\.\d+)?$` | `refseq_protein` |
| NCBI protein | `^[A-Z]{3}\d{5,7}(\.\d+)?$` or GenBank/GI | `ncbi_protein` |
| InChIKey | `^[A-Z]{14}-[A-Z]{10}-[A-Z]$` | `chemical` |
| InChI | starts with `InChI=` | `chemical` |
| SMILES | parses as SMILES; contains bond/atom tokens, not a gene-like word | `chemical` |
| PubChem CID | bare integer (context: chemical field) | `chemical` |
| Variant | `^[A-Za-z0-9]+ [A-Z]\d+[A-Z]$` (e.g. `BRAF V600E`) or `acc:p.Xaa#Xaa` | `variant` |
| Gene symbol | short all-caps alnum token matching HGNC-like shape | `gene_symbol` |
| Protein name / organism+name | free text | `protein_name` |
| List | newline/comma-separated multiple of the above | `batch` |
| Pair | protein-ish + chemical-ish separated by `+`/`,`/whitespace with a chem token | `protein_chemical_pair` |

`interpret(query) -> {query_type, tokens, confidence, needs_disambiguation}`.
When two rules plausibly match (e.g. `MDM2` is a gene *and* a name), both
candidate interpretations are returned ranked, and Stage-2 lookups decide.

## 5.2 Candidate disambiguation

For `gene_symbol`, `protein_name`, `organism+name`, or an ambiguous accession,
resolve to **candidate proteins** and present:

`{ accession, organism, taxon_id, gene, protein_name, length, review_status, has_structure }`

- Gene/name → UniProtKB search
  `?query=gene:BRAF AND organism_id:9606&fields=accession,id,protein_name,length,reviewed,organism_name`.
- Reviewed (Swiss-Prot) candidates rank above unreviewed (TrEMBL).
- Human (or a user-set organism filter) ranks first, but other organisms are
  still shown — never dropped.
- The user (or an explicit `taxon`/`accession` param) picks; the platform does
  **not** auto-collapse.

## 5.3 Cross-reference graph (Stage 2)

Build a directed multigraph of identifiers for the chosen entity. Nodes are
`(namespace, id)`; edges are `same_entity`, `isoform_of`, `sequence_identical`,
`maps_to`, each carrying its source.

Anchor lookups per entry point:

| Start | Bridge lookups |
|---|---|
| UniProt acc | UniProt entry → `uniParcId`, `uniProtKBCrossReferences` (RefSeq, PDB, Ensembl, GeneID, AlphaFold implied by acc) |
| UniParc UPI | UniParc entry → `uniParcCrossReferences[]` back to UniProt/RefSeq/EMBL (active + obsolete) |
| RefSeq / NCBI | E-utils `efetch` FASTA → CRC64 → UniParc `checksum:` search → UniProt; `elink` → Gene |
| PDB id | `rcsb.fetch_uniprot_accessions` (existing) → UniProt → rest |
| Raw FASTA | compute CRC64/MD5 → UniParc `checksum:` search → UniProt/RefSeq; if no exact hit, sequence-only record with a `novel_sequence` flag |
| Gene symbol | UniProt search (5.2) → chosen acc → rest |
| PubChem chem | independent compound graph (CID ↔ InChIKey ↔ parent/salt) |

### Sequence checksum as the identity key

CRC64 (UniProt/UniParc) and MD5 are computed for the query sequence and every
retrieved sequence. **Two records are the same *sequence*** iff checksums match;
they are the **same *protein* record** iff *checksum + taxon + isoform role*
match. This is what lets a raw FASTA resolve to a UniProt entry, and what stops a
human and mouse ortholog (same-ish sequence, different taxon) from merging.

## 5.4 Identity rules (what may and may **not** merge)

**MUST NOT merge on** gene symbol or protein name alone.

**Merge into one ProteinRecord only when all hold:**
1. Sequence checksum equal (or one is a documented isoform of the other), **and**
2. `taxon_id` equal, **and**
3. isoform role compatible (canonical↔canonical, or isoform explicitly linked via
   UniProt ALTERNATIVE PRODUCTS).

**Always kept distinct (separate records, linked by edges):**
- Isoforms (`isoform_of`), paralogs (`paralog_of`), orthologs (`ortholog_of`).
- PDB chains/entities (linked via `has_structure`, `chain_mapping`), and protein
  fragments/constructs (linked, with coverage).
- HLA/MHC alleles — **never** collapsed under one gene symbol; allele-level id
  (`HLA-A*02:01`) is preserved as the record identity (immunology workspace).

**Conflicts are reported, not resolved:** if UniProt and RefSeq disagree on
sequence, both are stored with a `sequence_conflicts` entry and provenance; the
UI shows the conflict rather than picking a winner (NFR-5).

## 5.5 Resolution output

```jsonc
{
  "query_type": "gene_symbol",
  "chosen": { "accession": "P15056", "taxon_id": 9606 },
  "candidates": [ /* all shown when ambiguous */ ],
  "xref_graph": {
    "nodes": [ {"ns":"UniProt","id":"P15056"}, {"ns":"UniParc","id":"UPI..."}, {"ns":"RefSeq","id":"NP_004324.2"}, {"ns":"PDB","id":"4MNE","qualifier":"A"}, {"ns":"AlphaFold","id":"AF-P15056-F1"} ],
    "edges": [ {"from":"UniProt:P15056","to":"UniParc:UPI...","rel":"sequence_identical","source":"UniParc"}, ... ]
  },
  "identity_key": { "taxon_id": 9606, "crc64": "...", "isoform": "canonical" },
  "warnings": [ "multiple organisms matched 'BRAF'; showing Homo sapiens (reviewed) — 12 others available" ]
}
```

## 5.6 Failure & edge behavior

- **No exact sequence match** for a raw FASTA → sequence-only record, `warnings:
  ["novel/unmatched sequence — no UniProt/UniParc record"]`, still fully
  analyzable (Stage 3).
- **Obsolete/replaced accession** → UniParc reveals the replacement; both the
  requested and current accession are recorded, requested one flagged obsolete.
- **Demerged accession** (one old acc → several new) → return all as candidates.
- **Rate-limited source** → partial graph with the missing bridge flagged
  `unavailable`, resolution still returns what it has.
