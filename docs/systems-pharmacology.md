# From single-site docking to systems-level questions

This note records where SnaCleX now sits on the path to multi-molecule,
multi-enzyme analysis, what the new panel engine can and cannot answer, and
what each further step actually requires. It exists so the gap between the
ambition and the current method is written down rather than implied.

## What shipped

`snaclex/panel.py` turns the workbench from a one-structure/one-site tool into
an **M ligands x N targets** engine:

- every ligand is docked into every target site in one run
  (`POST /api/jobs` with `kind: "panel"`), with live progress;
- one scoring grid per target, reused across all ligands, so cost is
  O(M*N) docks but only N grid builds;
- **per-target normalization**, because raw grid scores are not comparable
  across targets (deeper pockets produce more negative scores for *any*
  ligand). Ligand efficiency is standardized within each target column
  (`z_target`, higher = better), and only those standardized values are
  compared across a ligand's row to give `selectivity_gap`;
- each cell optionally carries **measured ChEMBL activity** for that
  (ligand, target) pair, so prediction sits next to experiment;
- one failed structure or compound degrades to an error cell, never an
  aborted run;
- `snaclex/systems.py` + `snaclex/data/systems/*.json` hold curated, versioned
  systems so a panel runs from an id. A catecholamine system ships as the
  first curation.

## What a panel can and cannot answer

**Can**: "Across this set of proteins and this set of molecules, which pairs
are geometrically and chemically plausible, which ligands look promiscuous,
and which look selective — ranked, reproducible, and carrying their method."
That is a triage instrument: it tells you what to measure next.

**Cannot**: predict a physiological event. A panel has no concentration, no
time, no expression level, no competition between ligands for the same site,
no flux, and no tissue. The scoring function is not calibrated to kcal/mol and
is not an affinity. Asking a panel to "predict a catecholamine storm" is a
category error in the same way asking a map for the traffic is: the map is a
real input to the question, and it is not the answer.

Three limits are worth stating plainly because they bite hardest on exactly
the pathway chosen as the first curation:

1. **Cofactors are stripped from the scoring grid.** `docking.build_grid`
   builds from `structure.protein_atoms` — standard amino acids only. For
   MAO-A/B the FAD is gone; for COMT the SAM and Mg are gone; for AADC the
   PLP is gone; for tyrosine hydroxylase the Fe and tetrahydrobiopterin are
   gone. Those sites are scored in an unphysical empty state. **This is the
   single largest correctness gap for enzyme panels** and is Rung 1 below.
2. **No protonation or tautomer handling.** Catecholamines are protonated
   amines at physiological pH. Docking their neutral heavy-atom geometry
   misstates precisely the electrostatics that drive binding in an anionic
   aromatic cage.
3. **Rigid ligand, rigid receptor.** One PubChem conformer, no induced fit.

## The capability ladder

Each rung is a real body of work with a different character. Rungs are ordered
by dependency: skipping ahead produces numbers that look like answers.

### Rung 1 — make a single cell trustworthy (prerequisite for everything)

A matrix of unreliable cells is an unreliable matrix, and normalization does
not fix it; it only makes the unreliability look orderly. Current measured
performance is 39–50% of redocked poses within 2 Å on a task the benchmark
README correctly notes is *easier* than the standard PoseBusters protocol
(crystal conformer retained, search centred on the true site). Work needed:

- **Retain cofactors in the grid.** Add an opt-in set of HETATM components to
  include as part of the receptor. Mechanically small, scientifically large.
- **Protonation at a stated pH**, at least for amines and carboxylates.
- **Ligand torsional flexibility**, or multiple PubChem conformers per ligand
  scored independently.
- **The optional Vina/GNINA track** already in `ROADMAP.md` Phase 4 — a panel
  is far more defensible scored by an engine reviewers already trust, with the
  dependency-free docker as the zero-install default.

### Rung 2 — calibrate, then report correlation not just rank

Panel z-scores are internal contrasts. To claim anything beyond "worth
measuring", each target needs a measured reference set (ChEMBL/BindingDB) and
each panel should report **Spearman correlation against measured affinity for
the subset of cells where a measurement exists**. The curated systems already
include positive controls for exactly this purpose: entacapone should top the
COMT column and selegiline the MAO-B column. If they do not, that panel's
rankings are unsupported and the result should say so rather than be read.

### Rung 3 — add the reaction network

A pathway layer — which enzyme acts on which metabolite, in which direction —
lets a panel's output be read structurally: "this metabolite engages the enzyme
two steps downstream of where it is produced". This is a graph over the
existing curation and needs no new physics. It is the natural next feature
after Rung 1 and the furthest a *structural* workbench should go on its own.

### Rung 4 — kinetics and dynamics

Turning "which pairs bind" into "what happens over time" requires Km, kcat,
Vmax, inhibition constants, enzyme concentrations and compartment volumes,
integrated as an ODE system. Those parameters are **measured** (BRENDA,
SABIO-RK), not derived from structure. Docking can at best supply priors on
whether an interaction exists at all.

This is where "predict the dynamics of a catecholamine surge" becomes a
well-posed question — and where it stops being a docking problem. The right
move is to **export to SBML and interoperate** with COPASI, Tellurium or
BioSimulators rather than reimplement a solver. SnaCleX's contribution is the
structural evidence feeding parameter choices, clearly labelled as such.

### Rung 5 — physiology and PK/PD

Tissue compartments, transport, clearance, receptor occupancy to functional
response. Out of scope for a structural workbench in any form; interface only.

**Recommended lane: Rungs 1–3.** They are achievable, they are what a
structural workbench is uniquely positioned to do, and Rung 1 is the one that
makes everything above it mean anything.

## Validation strategy

The panel feature should not be published on the strength of the panel feature.
The claim that needs evidence is *per-cell accuracy*, and the existing
benchmark harness already measures it. Concretely:

1. Finish Rung 1, then re-run `benchmarks/` and report the delta. Cofactor
   retention in particular should be reported as a before/after on the subset
   of complexes with cofactors.
2. Report against an external baseline (fpocket/P2Rank for pockets, Vina/smina
   for poses) on the same Astex/PoseBusters sets. A reviewer reads 39–50% as
   unanchored without one.
3. For panels specifically, report positive-control recovery: across curated
   systems, how often does the known inhibitor top its target's column?
4. Verify every curated system (`python -m snaclex.systems verify <id> --write`)
   and ship the verification evidence in the repository. Unverifiable
   assertions in a manuscript are what triggered the last screening decline;
   a committed verification record converts them into checkable ones.

## Where ChemTool fits

ChemTool's stated goal — predicting spectra and instrument readings for
designed molecules and metabolites — is a **cheminformatics problem, not a
structural one**, and the two projects should stay separate codebases. Their
natural relationship is a loop:

```
SnaCleX                              ChemTool
  structure-level hypothesis   -->   predicted MS/MS, NMR, retention
  (which molecule, which site)
             ^                                    |
             |                                    v
  re-prioritize from evidence  <--   match against measured metabolomics
```

The one thing worth deciding early is the **exchange format**, because getting
it wrong couples the projects permanently. Suggested: a small JSON record keyed
on **InChIKey** (not name, not CID — names are ambiguous and CIDs are
registry-specific), carrying SMILES, formula, charge, and a provenance block;
each tool attaches its own results under a namespaced key. That keeps both
sides independently useful and lets either be replaced.

Two cautions from the structural side, since they will apply equally to
ChemTool:

- **Predicted spectra are predictions.** The same provenance discipline
  SnaCleX applies to docking scores (method, parameters, version, measured
  comparison, limitations) should be built into ChemTool from the first
  commit, not retrofitted. It is far cheaper now.
- **Validate against a held-out measured set early** (e.g. MassBank, GNPS, or
  an in-house library) and report accuracy as a headline number. A spectral
  predictor with no reported accuracy is not usable by anyone else, and that
  is the gap that makes work hard to publish.

---

*Research-only. Nothing in SnaCleX, at any rung described here, constitutes
clinical guidance.*
