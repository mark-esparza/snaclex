# 03 — Source-by-Source Integration Table

One adapter per source. Endpoints below are the **documented public** endpoints.
Fields marked **`⚠ VERIFY`** were not round-tripped against a live response in
this offline environment and must be confirmed before production reliance (per
assumption **A1**). "Live-fetch only" means no bulk redistribution is assumed.

## Primary sources

### 1. RCSB PDB — experimental structures  ✅ *existing adapter `rcsb.py`*

| Aspect | Detail |
|---|---|
| Base | Data API `https://data.rcsb.org/rest/v1/core/...`; GraphQL `https://data.rcsb.org/graphql`; Search `https://search.rcsb.org/rcsbsearch/v2/query`; files `https://files.rcsb.org/download/{ID}.pdb|.cif` |
| Used for | Entry metadata, method, resolution; polymer entities/chains; bound ligands & chem components; UniProt xref per entity; full-text search; structural-similarity & sequence search; ModelServer for assemblies/ligand environments |
| Key fields | `struct.title`, `exptl[].method`, `rcsb_entry_info.resolution_combined`, `polymer_entities[].rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers[{database_name,database_accession}]`, `nonpolymer_entities`, auth/label residue numbering |
| Structure-by-UniProt | Search API filter on `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession` = UniProt acc, `database_name` = "UniProt". **`⚠ VERIFY`** exact attribute path for the current schema. |
| Sequence↔structure map | SIFTS via `rcsb_polymer_entity_align` / Sequence Coordinates API `https://sequence-coordinates.rcsb.org/graphql`. **`⚠ VERIFY`** |
| Auth / limits | Anonymous; be polite (batch via GraphQL). No key. |
| License | PDB data are public-domain (CC0-like); attribute wwPDB. |

### 2. UniProtKB — primary annotation  🆕 `uniprot.py`

| Aspect | Detail |
|---|---|
| Base | `https://rest.uniprot.org/uniprotkb/{accession}.json`; search `https://rest.uniprot.org/uniprotkb/search?query=...&format=json&fields=...` |
| Used for | Canonical sequence, isoforms, review status, protein/gene names, organism+taxId, function, catalytic activity, cofactors, localization, tissue, domains/regions, active/binding sites, signal peptide, TM, PTMs, natural variants, disease, interactions, xrefs, evidence codes |
| Key fields | `primaryAccession`, `uniProtkbId`, `entryType` (`"UniProtKB reviewed (Swiss-Prot)"` \| `"UniProtKB unreviewed (TrEMBL)"`), `sequence.{value,length,molWeight,crc64,md5}`, `organism.{scientificName,taxonId,lineage}`, `genes[].geneName.value` + `synonyms`, `proteinDescription.recommendedName.fullName.value`, `comments[]` (typed: FUNCTION, CATALYTIC ACTIVITY, COFACTOR, SUBCELLULAR LOCATION, TISSUE SPECIFICITY, ALTERNATIVE PRODUCTS→isoforms, DISEASE, INTERACTION), `features[]` (`type`,`location.{start,end}.value`,`description`,`featureCrossReferences`), `uniProtKBCrossReferences[]` (`database`,`id`,`properties`) |
| Isoforms | `comments[type=ALTERNATIVE PRODUCTS].isoforms[]` gives isoform IDs (`P38398-2`); isoform sequence via `.../uniprotkb/{iso}.fasta`. **`⚠ VERIFY`** isoform JSON shape. |
| Reviewed vs unreviewed | `entryType` string above — always surfaced (Swiss-Prot vs TrEMBL). |
| Auth / limits | Anonymous; supports pagination cursors; polite rate. No key required. |
| License | CC BY 4.0 — attribute UniProt Consortium. |

### 3. UniParc — sequence archive / provenance  🆕 `uniparc.py`

| Aspect | Detail |
|---|---|
| Base | `https://rest.uniprot.org/uniparc/{upi}.json`; search `https://rest.uniprot.org/uniparc/search?query=...` (e.g. `checksum:{CRC64}`, `uniprot:{acc}`, `upi:{UPI}`) |
| Used for | Sequence-level dedup, stable sequence identity (UPI), tracking identical sequences across DBs, historical versions, obsolete/replaced accessions |
| Key fields | `uniParcId` (`UPI0000...`), `sequence.{value,length,crc64,md5}`, `uniParcCrossReferences[]` (`database`,`id`,`active`,`versionI`,`version`,`created`,`lastUpdated`) |
| Role boundary | Sequence archive + provenance **only** — never a functional-annotation source (per task). |
| Auth / limits | Anonymous. License CC BY 4.0. |

### 4. NCBI Protein / RefSeq  🆕 `ncbi.py`

| Aspect | Detail |
|---|---|
| Base | E-utilities `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{esearch,esummary,efetch,elink}.fcgi`; optional NCBI Datasets v2 `https://api.ncbi.nlm.nih.gov/datasets/v2alpha/...` |
| Used for | Protein sequence records, reference proteins, genome-linked records, gene→protein, alternative products, organism/strain, CDS provenance, links to Gene/nucleotide/taxonomy/PubMed/CDD/variation |
| Key ops | `efetch?db=protein&id={acc}&rettype=fasta&retmode=text` (sequence); `esummary?db=protein&id={acc}&retmode=json` (title, organism, taxid, slen); `esearch?db=protein&term=...`; `elink?dbfrom=protein&db=gene\|nuccore\|pubmed&id=...` |
| RefSeq acc patterns | `NP_` (curated protein), `XP_` (predicted), `WP_` (non-redundant bacterial), `YP_`, `AP_`; versioned `NP_000537.3` |
| Auth / limits | 3 req/s anonymous, 10 req/s with `api_key` (`NCBI_API_KEY` env). Respect per-source limiter. |
| License | Public-domain US-gov data; cite NCBI. E-utilities usage policy applies. |

### 5. PubChem  ✅ *existing adapter `pubchem.py`* (extended for BioAssay)

| Aspect | Detail |
|---|---|
| Base | PUG-REST `https://pubchem.ncbi.nlm.nih.gov/rest/pug/...`; PUG-View `https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON` |
| Used for | CID/SID/name/synonym resolution; canonical & isomeric SMILES; InChI/InChIKey; formula/MW/charge; HBD/HBA; rotatable bonds; TPSA; XLogP; 2D/3D conformers; parent/salt/mixture/stereo relations; **BioAssay records, targets, activity outcomes, dose-response** |
| Key ops (existing) | `/compound/{ns}/{id}/property/{list}/JSON`; `/compound/cid/{cid}/SDF?record_type=3d\|2d`; PNG image |
| Key ops (new) | InChIKey `/compound/cid/{cid}/property/InChIKey,InChI/JSON`; parent `/compound/cid/{cid}/cids/JSON?cids_type=parent`; assay summary `/compound/cid/{cid}/assaysummary/JSON`; targets `/compound/cid/{cid}/assaysummary/JSON` → `Target GI/GeneID`; assay detail `/assay/aid/{aid}/description/JSON`. **`⚠ VERIFY`** assaysummary column names for activity value/units/outcome. |
| Auth / limits | No key; ≤5 req/s, avoid >400 req/min (PUG throttling headers). Existing `RateLimitError` handling applies. |
| License | Public-domain aggregate; individual assay depositors retain terms — cite AID + depositor. |

### 6. AlphaFold DB — predicted structures  🆕 `alphafold.py`

| Aspect | Detail |
|---|---|
| Base | API `https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}`; files `https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.{pdb,cif,bcif}` |
| Used for | Model id, UniProt mapping, predicted structure, per-residue pLDDT (B-factor column), PAE (via `paeImageUrl`/`paeDocUrl`), low-confidence/disordered regions, coverage, fragmentation of long proteins, model version/source |
| Key fields | `entryId`, `uniprotAccession`, `uniprotStart`, `uniprotEnd`, `pdbUrl`, `cifUrl`, `bcifUrl`, `paeImageUrl`, `paeDocUrl`, `modelCreatedDate`, `latestVersion`, `allVersions`, `globalMetricValue` (mean pLDDT). **`⚠ VERIFY`** `globalMetricValue` semantics and PAE JSON availability. |
| Confidence | pLDDT bands: ≥90 very high · 70–90 confident · 50–70 low · <50 very low. Never dock against low-confidence regions without prep assessment (Stage 4). |
| Guardrail | Never presented as experimental. Always badged "Predicted (AlphaFold vN)". |
| Auth / limits | Anonymous; be polite. License CC BY 4.0 (AlphaFold DB / EMBL-EBI + DeepMind). |

## Recommended enrichment sources (optional, independent adapters — Phase 2/3)

Each is **opt-in and env-gated**; failure or licensing never breaks the core.
"Verify before implementing" per the task.

| Source | Base | Used for | License / note (**`⚠ VERIFY`**) |
|---|---|---|---|
| **InterPro / InterProScan** | `https://www.ebi.ac.uk/interpro/api/` | Families, domains, repeats, sites | CC BY 4.0; InterProScan is heavy/opt-in |
| **Reactome** | `https://reactome.org/ContentService/` | Pathways, reactions | CC BY 4.0 |
| **IntAct** | `https://www.ebi.ac.uk/intact/ws/` (PSICQUIC) | Curated molecular interactions | CC BY 4.0 |
| **Gene Ontology / QuickGO** | `https://www.ebi.ac.uk/QuickGO/services/` | MF / BP / CC terms + evidence codes | CC BY 4.0 |
| **ClinVar** | E-utilities `db=clinvar` + VCV/RCV | Clinically interpreted variants + review status | Public; show review status, never interpret |
| **Ensembl / VEP** | `https://rest.ensembl.org/` | Genomic→protein consequence mapping | Apache-2.0 data terms; rate-limited |
| **STRING** | `https://string-db.org/api/` | PPI networks + confidence | CC BY 4.0; scores are predictions |
| **BioGRID** | `https://webservice.thebiogrid.org/` | Curated interactions | Requires access key |
| **IEDB** | `https://query-api.iedb.org/` | Epitopes (immunology workspace) | CC BY 4.0 (**`⚠ VERIFY`** API) |
| **PubMed** | E-utilities `db=pubmed` | Literature evidence | Public metadata; abstracts per publisher terms |

## Cross-cutting adapter contract

Every adapter must:

1. Accept a normalized identifier, return **normalized fragments + provenance**
   (`{value, source, source_id, source_version, retrieved_utc, evidence_category,
   limitations}`), never raw source JSON to callers.
2. Route all HTTP through the shared client (retry/backoff, conditional GET,
   per-source rate limit, API-key injection, cache).
3. Degrade gracefully: a miss/timeout/license-block returns "unavailable", not an
   exception that fails the whole record.
4. Be independently unit-testable via pure parse helpers over recorded fixtures
   (offline), matching the existing `test_*` convention.
5. Never fabricate a field. If a field is absent upstream, it is `null` with a
   provenance note — not inferred.
