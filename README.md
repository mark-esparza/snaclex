# SnaCleX

**Single-protein structure & ligand interaction analyzer — down to the atom.**

A research-only, dependency-free web tool. Load any experimental structure from
the **RCSB Protein Data Bank**, pick a bound molecule (drug, chemical, ion, or
metal), and SnaCleX computes and visualizes exactly how that molecule contacts
the protein at the atomic level — hydrogen bonds, hydrophobic packing, salt
bridges, metal coordination, and aromatic contacts — then cross-references the
chemical against **PubChem** and **ChEMBL**.

> Research-only. Interactions are geometric heuristics computed from experimental
> coordinates (no explicit hydrogens, no energy minimization). Not for clinical use.

> **Copyright © 2026 Mark Esparza. All rights reserved.** SnaCleX is proprietary
> software — see [LICENSE](LICENSE). No license to use, copy, modify, or
> distribute is granted; the source is public for reference only. Third-party
> components and data sources retain their own terms — see [NOTICE](NOTICE).

## What it does

1. **Load a protein** by PDB ID (or full-text search the PDB).
2. **Pick a bound molecule** — ligands/ions/metals are auto-detected from the structure.
3. **Atomic interaction profile** — every heavy-atom contact between the molecule
   and the protein is classified by geometry and listed with distances, plus a
   per-residue "binding hot spot" summary.
4. **Pocket detection** — find geometric cavities (LIGSITE) so you can target
   *apo* structures with no bound ligand; pockets are ranked by a druggability
   score and listed with volume, enclosure, and lining residues (see Pockets below).
5. **Dock a PubChem chemical into a site** — place a real 3D conformer of any
   drug/chemical into a bound-ligand site *or a detected pocket* via Monte-Carlo
   rigid-body search, then profile the predicted pose's atomic interactions.
6. **3D atomic viewer** (3Dmol.js) — protein cartoon + ligand sticks + dashed
   interaction lines, optional molecular surface; docked poses render as black
   sticks and detected pockets as a translucent sphere over their lining residues.
7. **Chemical lookup** (PubChem) — formula, MW, XLogP, TPSA, H-bond donors/acceptors,
   rotatable bonds, SMILES, Lipinski druglikeness, plus ChEMBL development status.
8. **Research report** — plain-language summary + generated hypotheses, exportable as `.txt`.

## Pockets (apo-structure docking)

When a structure has no bound ligand, SnaCleX finds cavities with a pure-Python
**LIGSITE** implementation: the protein is placed on a grid, grid points buried
inside the molecular volume are flagged by scanning 7 directions for
protein-solvent-protein enclosure, and the buried points are clustered into
ranked cavities. Each pocket reports a centroid, volume, mean enclosure (5–7),
lining residues, and a 0–100 druggability score (volume + enclosure + lining).
You can **Dock here** straight from any detected pocket.

Validation: detecting pockets in 1HSG/1CA2 (using protein atoms only) recovers
the true ligand site as the top-ranked pocket (centroid within ~1–5 Å of the
crystallographic ligand).

## Docking (true pose search)

AutoDock Vina/RDKit won't install cleanly on Windows + Python 3.14, so SnaCleX
ships an AutoDock-style **grid-map docker** written in pure Python:

- The ligand's real 3D conformer comes from **PubChem** (`record_type=3d`), so no
  conformer-generation dependency is needed.
- A scoring grid is precomputed over the pocket with **steric**, **hydrogen-bond**,
  and **hydrophobic** channels (trilinear-interpolated AutoDock-style affinity maps).
- The rigid ligand is searched with **Monte-Carlo translation/rotation** plus
  simulated-annealing acceptance; the best pose is returned and its atomic
  interactions are profiled with the same engine used for crystal ligands.

It is approximate and **research-only**: rigid ligand, empirical score (lower =
better fit, *not* kcal/mol), no explicit solvent or full force field. Validation:
redocking benzamidine into trypsin (`3PTB`) reproduces the known S1-pocket pose —
top contact **Asp189** with salt bridges, at **~2.0 Å redock RMSD** (sub-2 Å is the
usual docking-success threshold), in ~2 s.

### Redock RMSD
When you dock a chemical into the site of an existing crystallographic ligand of
the same size, SnaCleX reports a **nearest-atom RMSD** between the predicted and
experimental pose — an automatic accuracy readout for redocking validation.

### Batch virtual screening
The Docking tab can screen a list of chemicals (up to 10) into one target. The
scoring grid is built **once** and reused for every ligand, so each additional
compound costs only the Monte-Carlo search (~2 s). Results are ranked by predicted
fit, with both total score and per-atom score (ligand efficiency) shown — total
score favours larger molecules, so per-atom is the fairer cross-size comparison.
`GET /api/screen?pdb=ID&chems=a,b,c&comp=INDEX` (or `&pocket=INDEX`).

## Interaction criteria (heavy-atom, no explicit H)

| Type | Geometric rule |
| --- | --- |
| Hydrogen bond / polar | N/O .. N/O ≤ 3.6 Å (S up to 3.9 Å) |
| Salt bridge | charged side chain (Asp/Glu/Arg/Lys/His) .. opposite-charge ligand atom ≤ 4.0 Å |
| Hydrophobic | C .. C, 2.8–4.0 Å (closest per residue) |
| Metal coordination | metal .. O/N/S ≤ 2.9 Å |
| Aromatic (possible π) | aromatic ring centroid .. nearest ligand heavy atom ≤ 4.5 Å |

These mirror a simplified PLIP-style approach so the tool runs with zero
scientific dependencies. They are screening heuristics, not force-field energies.

## Run

```bash
python server.py            # http://127.0.0.1:8010
python server.py --port 8000
```

Pure Python standard library — no pip installs required. The frontend loads
3Dmol.js from a CDN. Internet access is needed for RCSB / PubChem / ChEMBL.

## Deploy to the cloud (free)

The app is a single stdlib Python process with **no dependencies**, so most PaaS
hosts run it as-is. It reads `$PORT` and binds `0.0.0.0` automatically when that
variable is set (see `Procfile`).

**Render (recommended, free tier):**
1. Push this repo to GitHub.
2. On [render.com](https://render.com) → *New* → *Web Service* → connect the repo.
3. Settings: **Build Command** `python --version` · **Start Command** `python server.py`.
4. Deploy. Render gives you a public `https://<name>.onrender.com` URL.

**Railway / Fly.io / Heroku** work the same way via the included `Procfile`.

The host needs outbound internet (for RCSB/PubChem/ChEMBL) — all the major
platforms allow this by default. Note: free tiers sleep when idle, so the first
request after a pause may take ~30 s to wake.

## Try it

- `1HSG` — HIV-1 protease + indinavir (MK1): rich H-bond + hydrophobic pocket; catalytic Asp25 is the top contact.
- `1CA2` — carbonic anhydrase II: zinc coordinated by His94/His96/His119.
- `3PTB` — trypsin + benzamidine + calcium.
- Chemical box: `aspirin`, `imatinib`, `zinc`.

## Layout

```
server.py              stdlib HTTP server + JSON API
snaclex/
  http_util.py         urllib helpers
  rcsb.py              PDB structure + metadata + search
  pubchem.py           compound lookup + Lipinski
  chembl.py            optional drug/bioactivity cross-reference
  pdbparse.py          dependency-free PDB parser
  interactions.py      atomic interaction profiler
  pockets.py           LIGSITE geometric cavity finder
  docking.py           grid-map Monte-Carlo rigid-body docker
  report.py            summary + hypothesis generator
web/
  index.html style.css app.js
legacy/                previous StructInteract CLI (archived)
```

## API

The full, live contract is served at `GET /api/docs` (rendered at `/api.html`).

| Endpoint | Returns |
| --- | --- |
| `GET /api/analyze?pdb=ID` | metadata, chains, components, raw PDB text |
| `GET /api/interactions?pdb=ID&comp=INDEX` | interaction profile + report |
| `GET /api/chemical?q=NAME` | PubChem properties + druglikeness + ChEMBL |
| `GET /api/pockets?pdb=ID` | detected geometric cavities (ranked) + methods/provenance |
| `GET /api/evolution?pdb=ID` | Pfam conservation per residue/pocket + methods/provenance |
| `GET /api/search?q=TEXT` | PDB full-text search results |
| `GET /api/version` · `GET /api/docs` | version · machine-readable API contract |
| `POST /api/jobs` → `GET /api/jobs/{id}` | submit a docking/screening job (`kind` = `dock`/`screen`) and poll its status/result |
| `POST /api/upload` | analyze a user-supplied PDB-format file (returns an upload id usable as `ID` above) |

`ID` is a 4-char PDB id **or** an upload id from `POST /api/upload`.

## Development

SnaCleX has a stdlib-only test suite (no pytest, no third-party deps). Run it
from the repo root:

```
python -m unittest discover -s tests -v
```

The tests are fully offline — the pure-compute modules run on small synthetic
structures, and the HTTP layer is exercised through a mock seam, so no live RCSB
/ PubChem / ChEMBL calls are made. CI (`.github/workflows/ci.yml`) byte-compiles
the sources and runs the suite on Python 3.11–3.13.

See [ROADMAP.md](ROADMAP.md) for planned work and [CHANGELOG.md](CHANGELOG.md)
for recent changes.

### Benchmarking the docker

A dependency-free redocking benchmark self-docks each crystallographic ligand
back into its own receptor and reports pose-recovery metrics (top-1 RMSD ≤ 2 Å,
median/mean RMSD):

```
python -m snaclex.benchmark path/to/*.pdb --out benchmark_results.json
# or fetch from RCSB by ID:
python -m snaclex.benchmark 1HSG 1CA2 --out benchmark_results.json
```

Point it at PoseBusters / CrossDocked / PDBbind structures to benchmark at
scale. If `benchmark_results.json` is committed at the repo root, the summary is
surfaced in every docking result's methods block. RMSD is a nearest-atom proxy
(see `docking.rmsd_to_reference`), so read sub-2 Å as "pose recovered".
