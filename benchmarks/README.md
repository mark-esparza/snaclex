# SnaCleX validation benchmark

Scripts and results behind the Validation section of the SnaCleX preprint
(runs of 5 October 2026, SnaCleX 0.1.0).

## Data

The Astex Diverse set (85 complexes) and the PoseBusters Benchmark set
(428 complexes), as prepared by Buttenschoen, Morris and Deane:
Zenodo record 8278563, file `posebusters_paper_data.zip`. Unzip it so that the
`astex_diverse_set/`, `posebusters_benchmark_set/` and `*_ids.txt` files sit in
a folder named `pb` next to the folder that holds these scripts and the
`snaclex` package (the scripts read `../pb`).

## What each script does

| Script | Environment | Output |
|---|---|---|
| `run_posebusters.py [workers]` | Python 3.11+, SnaCleX only | `<set>_results.jsonl`: pocket recovery, contacts, and rigid redocking of the crystal ligand conformer under seeds 0, 1, 2; docked poses in `poses/` |
| `plip_compare.py [workers]` | conda env with `plip`, `openbabel` (see `run_plip.ps1`) | `plip_<set>.jsonl`: residue sets per interaction type from PLIP and SnaCleX |
| `analyze.py all` | same env plus `rdkit` | `summary.json` and `<set>_with_sym.json` (adds symmetry corrected RMSD via RDKit CalcRMS) |
| `timing.py` | Python 3.11+, SnaCleX only | `timing.json`: single process timings on 20 Astex complexes |
| `run_example.py`, `make_figs.py`, `shoot.py`, `fig_template.html` | SnaCleX; Microsoft Edge for headless screenshots | trypsin and benzamidine (3PTB) worked example and Figures 2 to 4 |
| `plot_rmsd_curves.py` | numpy, matplotlib | Figure 5 |

## Headline results

| | Astex Diverse (85) | PoseBusters (428) |
|---|---|---|
| Redocked within 2 Å, single run (mean of 3 seeds) | 39.2% | 40.0% |
| Redocked within 2 Å, best score of 3 runs | 50.6% | 47.2% |
| Pocket centroid within 6 Å of the ligand | 50.6% | 51.4% |
| PLIP contact residues recovered (any type) | 97.4% | 95.5% |
| Precision against PLIP (any type) | 61.3% | 64.6% |

RMSD is the symmetry corrected heavy atom RMSD without superposition. The
redocking task keeps the crystal conformer and centers the search on the true
site, so it is easier than the standard PoseBusters protocol. PLIP could not
finish 8F4J (about 58,000 protein atoms), so the PLIP comparison covers 427
PoseBusters complexes. Full numbers are in `results/summary.json`.
