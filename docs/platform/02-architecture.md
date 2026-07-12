# 02 — System Architecture & Data Flow

Two views of the same system: the **target production architecture** (recommended
stack, the destination at scale) and the **vertical-slice architecture** (what
ships now, stdlib-first). They share the same layer boundaries so the slice
implementation can be swapped for the production one behind stable interfaces.

## 2.1 System architecture (target)

```mermaid
flowchart TB
    subgraph Clients
        UI[React + TS web app<br/>Mol* viewer]
        CLI[Python client / CLI]
    end

    subgraph API[API layer — FastAPI + Pydantic, versioned /api/v1]
        GW[Request routing · auth · rate limit · OpenAPI]
    end

    subgraph Services[Domain services]
        IDR[Identifier Resolution]
        SEQ[Sequence Analysis]
        STR[Structural Analysis<br/>«existing SnaCleX modules»]
        CHEM[Chemical Normalization]
        EVID[Evidence Integration Engine]
        SEARCH[Search / comparison]
    end

    subgraph Adapters[Source-adapter layer «one per source, isolated»]
        UP[UniProtKB]
        UPI[UniParc]
        NCBI[NCBI Protein / RefSeq]
        RCSB[RCSB PDB]
        AF[AlphaFold DB]
        PC[PubChem]
        OPT[(optional: InterPro · Reactome · IntAct · GO ·<br/>ClinVar · Ensembl · STRING/BioGRID · PubMed · IEDB)]
    end

    subgraph Infra[Cross-cutting]
        HTTP[HTTP client: retry/backoff · conditional GET · per-source rate limit · API keys]
        CACHE[(Response cache<br/>disk / Redis)]
        JOBS[Job queue<br/>threads / Celery+Redis]
    end

    subgraph Storage
        PG[(PostgreSQL<br/>normalized records + evidence)]
        OBJ[(Object storage<br/>structures · alignments · reports)]
        IDX[(Search index)]
        GRAPH[(optional graph DB<br/>knowledge graph)]
    end

    UI --> GW
    CLI --> GW
    GW --> IDR & SEQ & STR & CHEM & EVID & SEARCH
    IDR & SEQ & STR & CHEM & EVID --> Adapters
    Adapters --> HTTP --> CACHE
    STR --> JOBS
    EVID --> JOBS
    Services --> PG
    STR --> OBJ
    SEARCH --> IDX
    EVID --> GRAPH
```

### Layer responsibilities

| Layer | Responsibility |
|---|---|
| **Source adapters** | One module per external source. Turn source-native JSON into normalized fragments with provenance. No cross-source logic. Independently testable, independently disableable. |
| **HTTP / infra** | Single choke point for retry, backoff, conditional requests, per-source rate limits, API-key injection, and caching. (Today: `snaclex/http_util.py` + `snaclex/cache.py`.) |
| **Identifier resolution** | Classify the query; build the cross-reference graph; decide protein identity; keep 1:N and N:N links. |
| **Sequence analysis** | Pure-compute sequence metrics + retrieved sequence features, merged into the record. |
| **Structural analysis** | The *existing* SnaCleX capability (`pdbparse`, `interactions`, `pockets`, `docking`, `evolution`) plus predicted-structure confidence handling. |
| **Chemical normalization** | Compound identity + PDB-ligand↔PubChem match-quality classification. |
| **Evidence integration** | Assemble the universal evidence objects, assign level A–F, rank, and forbid flattening. |
| **Storage** | Normalized records + evidence (relational), large binaries (object store), search index, optional graph. |
| **Jobs** | Docking, large alignments, batch, prediction — never inline in a request. |

## 2.2 Vertical-slice architecture (what ships now)

Same boundaries, stdlib implementations. No FastAPI, Postgres, Redis, or React;
everything runs from `python server.py` with an empty `requirements.txt`.

```mermaid
flowchart TB
    UI[web/ vanilla-JS SPA + 3Dmol.js]
    SRV[server.py — stdlib ThreadingHTTPServer + JSON API]

    subgraph New[New slice modules «snaclex/»]
        IDR[idresolve.py]
        REC[proteinrecord.py]
        SEQ[seqanalysis.py]
        EV[evidence.py]
    end
    subgraph Adapt[Adapters «snaclex/»]
        UP[uniprot.py]
        UPI[uniparc.py]
        NC[ncbi.py]
        AF[alphafold.py]
        RC[rcsb.py ✓existing]
        PC[pubchem.py ✓existing]
    end
    subgraph Struct[Structural Analysis «existing, preserved»]
        PP[pdbparse] --- IN[interactions] --- PK[pockets] --- DK[docking] --- EVv[evolution]
    end
    HU[http_util.py + cache.py] 
    JB[jobs.py]

    UI --> SRV --> IDR & REC & SEQ & EV
    IDR --> UP & UPI & NC & RC & AF
    REC --> UP & UPI & NC & AF
    SEQ --> UP
    EV --> RC & PC
    SRV --> Struct
    Adapt --> HU
    SRV --> JB
```

**Interface stability:** each slice module exposes a plain-function/plain-dict
API (no framework types). Moving to FastAPI later means wrapping these functions
in Pydantic response models and swapping `http_util`'s cache for Redis and
`jobs.py` for Celery — callers do not change.

## 2.3 Data-flow diagram — "resolve a query to a report"

```mermaid
sequenceDiagram
    autonumber
    actor R as Researcher
    participant API
    participant IDR as Identifier Resolution
    participant AD as Source Adapters
    participant SEQ as Sequence Analysis
    participant STR as Structure Availability
    participant EV as Evidence Engine
    participant CACHE

    R->>API: query (accession | gene | FASTA | chem | pair | variant)
    API->>IDR: interpret(query)
    IDR-->>API: query_type + candidates (organism, acc, len, review)
    Note over R,API: if ambiguous → return candidates, ask user to pick
    API->>IDR: resolve(chosen)
    IDR->>AD: UniProt / UniParc / NCBI / RCSB lookups
    AD->>CACHE: get-or-fetch (retry, backoff, cond. GET)
    CACHE-->>AD: cached bytes | live JSON
    AD-->>IDR: normalized fragments + provenance
    IDR-->>API: ProteinRecord (identity resolved, xref graph)
    API->>SEQ: analyze(record.sequence)
    SEQ-->>API: composition, MW, pI, charge@pH, GRAVY, low-complexity (all "calculated")
    API->>STR: assess(record)
    STR->>AD: RCSB by-UniProt search + AlphaFold prediction API
    STR-->>API: structure hierarchy (exp | homolog | predicted | none) + coverage/identity/confidence
    opt protein–chemical pair
        API->>EV: evidence(protein, chemical)
        EV->>AD: PDB bound-ligand · PubChem BioAssay · (curated · homology)
        EV-->>API: evidence[] typed A–F, ranked, never merged
    end
    API-->>R: report (13 tabs) + provenance + downloadable snapshot
```

## 2.4 Data-flow diagram — "ingestion vs interactive" (NFR-4)

Reports must not fan out live calls to every source on every open. Three
separate paths share the cache and record store:

```mermaid
flowchart LR
    subgraph Interactive[Interactive query path «low latency»]
        Q[query] --> RECget[read ProteinRecord]
        RECget -->|hit| RPT[render report]
        RECget -->|miss / stale| ING
    end
    subgraph Ingest[Ingestion path «on-demand or batch»]
        ING[resolve + fetch adapters] --> NORM[normalize + provenance] --> STORE[(record store + cache)]
        STORE --> RPT
    end
    subgraph Refresh[Cache-refresh path «scheduled»]
        CRON[release tracking / TTL] --> REFR[conditional re-fetch changed sources] --> STORE
    end
```

- **Interactive** reads a materialized record; only a cold/stale miss triggers
  ingestion.
- **Ingestion** does the fan-out once, normalizes, stamps provenance, persists.
- **Refresh** watches source releases (UniProt release, PDB weekly, PubChem)
  and re-fetches only changed records with conditional requests.
