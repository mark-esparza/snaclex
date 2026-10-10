"""Many-to-many ligand x target panels (systems-level interaction matrices).

Everything else in SnaCleX answers a *one structure, one site* question. A
physiological question -- "which of these circulating metabolites engage which
of these enzymes and receptors, and how selectively?" -- is a **matrix**
question: M ligands against N targets, scored on a common footing.

This module is that matrix engine. It is deliberately I/O-free: structure
loading, ligand fetching, site resolution and grid construction are all injected
by the caller (``server.py`` passes its cached loaders; tests pass fakes), so the
scoring and normalization logic is unit-testable with no network.

Why normalization is the whole point
------------------------------------
Raw docking scores from :mod:`snaclex.docking` are **not comparable across
targets**. The grid score is an unnormalized sum over ligand heavy atoms inside
a box whose occupancy differs per site, so a deep buried pocket produces more
negative numbers than a shallow surface groove for *any* ligand. Reading a raw
matrix row-wise would therefore rank targets by pocket burial, not by affinity.

So the engine standardizes **within each target column** before anything is
compared across targets:

1. ``ligand_efficiency`` (score / heavy atom) is the per-cell comparable
   quantity -- it removes the first-order dependence on ligand size, which spans
   a wide range across a metabolite panel.
2. ``z_target`` standardizes ligand efficiency within one target's column, sign
   flipped so **higher is better**. This answers "is this ligand unusually good
   *for this target*, relative to the rest of the panel?"
3. Only those already-standardized values are compared across a ligand's row, to
   give ``selectivity_gap`` -- best z_target minus runner-up. This answers "does
   this ligand prefer one target, or hit everything equally?"

Both readouts are *relative to the panel you submitted*. Adding or removing
ligands changes every z-score. They are panel-internal contrasts, never absolute
affinities.

What this is not
----------------
This does not simulate a physiological event. There is no concentration, no
time, no expression level, no competition between ligands for the same site, no
flux through a pathway. A panel is a static, structure-derived *hypothesis
matrix*: a ranked, provenance-carrying starting set of "this pair is worth
measuring", to be confirmed against binding data (see ``measured`` cells, which
carry ChEMBL activity when the caller supplies a lookup) and wet-lab assays.
See ``docs/systems-pharmacology.md`` for what would be required to go further.
"""

from __future__ import annotations

import datetime
import math
import statistics

from . import __version__, docking, interactions
from .pdbparse import Structure

# A panel holds every target's receptor in memory at once, so a 16-target run
# of large assemblies would be hundreds of MB. Nothing in a panel reads protein
# atoms far from the docking site: a docked atom cannot leave the grid box
# (GRID_HALF from the center), and the widest interaction cutoff is ~4.5 A on
# top of that. Keeping whole residues within this radius is therefore exact,
# not an approximation -- and it is a radius, not a cutoff, so the margin below
# is pure safety. Residues (not atoms) are kept whole so aromatic rings survive
# intact for ring-centroid detection.
SITE_RADIUS = docking.GRID_HALF + 12.0

# Guard the invariant above: the retained shell must comfortably exceed the
# furthest reach of any scored contact.
assert SITE_RADIUS > docking.GRID_HALF + max(
    interactions.HB_MAX, interactions.SALT_MAX, interactions.HYDRO_MAX,
    interactions.METAL_MAX, interactions.ARO_CENTROID_MAX,
)

# Bounds that keep a panel's runtime sane. A panel is O(M*N) docks plus N grid
# builds; at the default search budget a full 12x16 grid is already minutes.
MAX_TARGETS = 16
MAX_LIGANDS = 24
MAX_CELLS = 192

# Panels use a lighter search than a single focused dock (which uses 220 seeds),
# because a panel trades per-pose precision for breadth. Raise via `seeds=` when
# a specific pair matters.
PANEL_SEEDS = 120
PANEL_MC_STEPS = 40

# Contact residues retained per cell (full profiles across 192 cells are huge).
TOP_RESIDUES = 6


class PanelError(ValueError):
    """Panel could not be set up (bad spec, nothing resolvable)."""


def _mean_sd(values):
    """Return (mean, population stdev) or (mean, None) when sd is undefined."""
    if not values:
        return None, None
    mean = statistics.fmean(values)
    if len(values) < 2:
        return mean, None
    sd = statistics.pstdev(values)
    return mean, (sd if sd > 1e-9 else None)


def _round(value, digits=3):
    return None if value is None else round(value, digits)


def _normalize(cells, targets, ligands):
    """Attach z_target / rank_in_target / selectivity to scored cells, in place.

    Returns the per-target normalization record (mean, sd, n) so the result can
    state exactly what each z-score was computed against.
    """
    by_target: dict[str, list] = {t["id"]: [] for t in targets}
    for cell in cells:
        if cell.get("error") is None:
            by_target[cell["target_id"]].append(cell)

    norm = {}
    for target in targets:
        tid = target["id"]
        column = by_target[tid]
        les = [c["ligand_efficiency"] for c in column]
        mean, sd = _mean_sd(les)
        norm[tid] = {
            "n_scored": len(column),
            "mean_ligand_efficiency": _round(mean),
            "sd_ligand_efficiency": _round(sd),
            # A z-score needs >=2 ligands and real spread; say so when absent
            # rather than silently emitting zeros.
            "z_available": sd is not None,
        }
        # Lower (more negative) ligand efficiency is a better score, so negate
        # to make a positive z mean "better than this target's panel average".
        for cell in column:
            cell["z_target"] = (
                None if sd is None
                else _round(-(cell["ligand_efficiency"] - mean) / sd)
            )
        for rank, cell in enumerate(
            sorted(column, key=lambda c: c["ligand_efficiency"]), start=1
        ):
            cell["rank_in_target"] = rank
            cell["n_in_target"] = len(column)

    # Selectivity: compare a ligand's *already standardized* values across its
    # row. Comparing raw scores here would rank pocket burial, not preference.
    selectivity = {}
    for ligand in ligands:
        lid = ligand["id"]
        row = [
            c for c in cells
            if c["ligand_id"] == lid and c.get("error") is None
            and c.get("z_target") is not None
        ]
        if not row:
            selectivity[lid] = {
                "best_target_id": None, "selectivity_gap": None,
                "note": "no target in this panel produced a usable z-score",
            }
            continue
        row.sort(key=lambda c: c["z_target"], reverse=True)
        best = row[0]
        gap = (
            None if len(row) < 2
            else _round(best["z_target"] - row[1]["z_target"])
        )
        selectivity[lid] = {
            "best_target_id": best["target_id"],
            "best_z": best["z_target"],
            "runner_up_target_id": row[1]["target_id"] if len(row) > 1 else None,
            "selectivity_gap": gap,
            "n_targets_scored": len(row),
        }
    return norm, selectivity


def _profile_summary(profile):
    """Compress a full interaction profile to what a matrix cell needs."""
    residues = [
        {
            "res_id": r["res_id"],
            "res_name": r["res_name"],
            "res_seq": r["res_seq"],
            "chain": r["chain"],
            "total": r["total"],
            "types": r["types"],
            "min_distance": r["min_distance"],
        }
        for r in profile["contact_residues"][:TOP_RESIDUES]
    ]
    return {
        "counts": profile["counts"],
        "interaction_total": profile["interaction_total"],
        "contact_residue_count": profile["contact_residue_count"],
        "top_residues": residues,
    }


def site_subset(structure, center, radius=SITE_RADIUS):
    """Return a Structure holding only whole residues near ``center``.

    Keeps a panel's memory proportional to its sites rather than to the size of
    its assemblies. Exact for panel purposes: see :data:`SITE_RADIUS`. Residues
    are kept whole so that aromatic-ring perception is unaffected.
    """
    cx, cy, cz = center
    r2 = radius * radius

    keep = set()
    for a in structure.protein_atoms:
        dx, dy, dz = a.x - cx, a.y - cy, a.z - cz
        if dx * dx + dy * dy + dz * dz <= r2:
            keep.add((a.chain, a.res_seq, a.icode))

    protein = [
        a for a in structure.protein_atoms
        if (a.chain, a.res_seq, a.icode) in keep
    ]
    if not protein:
        # Nothing within reach: hand back the original rather than an empty
        # receptor, so the failure shows up as a bad score, not a silent void.
        return structure
    return Structure(
        atoms=protein,
        protein_atoms=protein,
        components=[],   # panels score docked poses, not deposited ligands
        chains=sorted({a.chain for a in protein}),
    )


def cofactor_atoms(structure, wanted):
    """Resolve declared cofactor residue names to components in `structure`.

    Returns ``(components, missing_names)``. Missing cofactors are reported
    rather than raised: a site scored without its cofactor is still scoreable,
    it is just less meaningful, and the caller surfaces that as a warning.
    """
    if not wanted:
        return [], []
    names = {str(n).upper() for n in wanted}
    found = [c for c in structure.components if c.res_name.upper() in names]
    missing = sorted(names - {c.res_name.upper() for c in found})
    return found, missing


def _check_box_fit(ligand_atoms):
    """Flag ligands too large for the fixed docking box (pose would be clipped)."""
    xs = [a["x"] for a in ligand_atoms]
    ys = [a["y"] for a in ligand_atoms]
    zs = [a["z"] for a in ligand_atoms]
    span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
    # The ligand centroid may wander TRANS_HALF from the site center, so the
    # furthest atom can reach TRANS_HALF + span/2 and must stay inside GRID_HALF.
    reach = docking.TRANS_HALF + span / 2.0
    return {
        "max_span_A": round(span, 2),
        "fits_box": reach <= docking.GRID_HALF,
        "note": (
            None if reach <= docking.GRID_HALF else
            f"ligand span {span:.1f} A may not fit the "
            f"{2 * docking.GRID_HALF:.0f} A box; pose is likely clipped"
        ),
    }


def run_panel(
    targets,
    ligands,
    *,
    load_structure,
    load_ligand,
    resolve_site,
    build_grid=None,
    measured_activity=None,
    progress=None,
    seeds=PANEL_SEEDS,
    mc_steps=PANEL_MC_STEPS,
    random_seed=0,
):
    """Dock every ligand into every target site and return a normalized matrix.

    Parameters
    ----------
    targets, ligands:
        Spec dicts. A target needs ``id`` and ``pdb`` plus a site hint consumed
        by ``resolve_site``; a ligand needs ``id`` and ``query``.
    load_structure(pdb) -> (text, structure, meta)
    load_ligand(query) -> {"atoms", "cid", "source", "name", "formula"}
    resolve_site(pdb, structure, target) -> (center, label)
    build_grid(pdb, structure, center, extra_atoms) -> Grid
        Optional; defaults to an uncached ``docking.build_grid``. The server
        passes its per-(structure, site, cofactors) cached builder.
    measured_activity(ligand, target) -> dict | None
        Optional lookup for experimental activity (e.g. ChEMBL) attached to each
        cell as ``measured``, so predictions sit next to real measurements.
    progress(done, total, label) -> None
        Optional callback invoked after each cell.

    One failing ligand or target degrades to an error cell / error row rather
    than aborting the panel -- a 200-cell run should not be lost to one 404.
    """
    if not targets:
        raise PanelError("A panel needs at least one target")
    if not ligands:
        raise PanelError("A panel needs at least one ligand")
    if len(targets) > MAX_TARGETS:
        raise PanelError(f"Too many targets (max {MAX_TARGETS})")
    if len(ligands) > MAX_LIGANDS:
        raise PanelError(f"Too many ligands (max {MAX_LIGANDS})")
    if len(targets) * len(ligands) > MAX_CELLS:
        raise PanelError(
            f"Panel too large: {len(targets)}x{len(ligands)} exceeds "
            f"{MAX_CELLS} cells"
        )
    if build_grid is None:
        def build_grid(_pdb, structure, center, extra_atoms=None):
            return docking.build_grid(structure, center, extra_atoms)

    total = len(targets) * len(ligands)
    done = 0

    # --- Resolve ligands once; a ligand is fetched one time for the whole panel.
    resolved_ligands = []
    for spec in ligands:
        entry = {
            "id": spec["id"],
            "query": spec.get("query") or spec["id"],
            "label": spec.get("label") or spec.get("query") or spec["id"],
            "role": spec.get("role"),
        }
        try:
            lig = load_ligand(entry["query"])
            entry.update({
                "cid": lig.get("cid"),
                "name": lig.get("name") or entry["label"],
                "formula": lig.get("formula"),
                "n_heavy_atoms": len(lig["atoms"]),
                "coord_source": lig.get("source"),
                "geometry": _check_box_fit(lig["atoms"]),
                "error": None,
            })
            entry["_atoms"] = lig["atoms"]
        except Exception as exc:  # noqa: BLE001 - degrade to an error row
            entry.update({"error": str(exc), "_atoms": None})
        resolved_ligands.append(entry)

    # --- Resolve targets once; the grid is the expensive part and is built per
    # target, not per cell. This is what makes M*N affordable.
    resolved_targets = []
    for spec in targets:
        entry = {
            "id": spec["id"],
            "pdb": spec.get("pdb"),
            "label": spec.get("label") or spec["id"],
            "role": spec.get("role"),
            "gene": spec.get("gene"),
        }
        try:
            _text, structure, meta = load_structure(entry["pdb"])
            center, site_label = resolve_site(entry["pdb"], structure, spec)
            entry.update({
                "title": meta.get("title"),
                "resolution_A": meta.get("resolution_A"),
                "experimental_method": meta.get("experimental_method"),
                "site": site_label,
                "center": [round(c, 2) for c in center],
                "error": None,
            })
            # Cofactors named by the curation join the rigid receptor;
            # without them a FAD/SAM/PLP-dependent site is scored as an empty
            # cavity (see docs/systems-pharmacology.md).
            cofactors, missing = cofactor_atoms(structure, spec.get("cofactors"))
            entry["cofactors_included"] = sorted(
                {c.res_name.upper() for c in cofactors}
            )
            entry["cofactors_missing"] = missing
            if missing:
                entry["warning"] = (
                    "declared cofactor(s) not found in this entry: "
                    + ", ".join(missing)
                    + " — this site is scored without them"
                )
            extra = [a for comp in cofactors for a in comp.atoms]

            # Build the grid from the full structure (the caller's cache keys
            # on the real entry), but retain only the site shell for the
            # per-cell interaction profiling that follows.
            entry["_grid"] = build_grid(entry["pdb"], structure, center, extra)
            entry["_structure"] = site_subset(structure, center)
            entry["n_site_residue_atoms"] = len(entry["_structure"].protein_atoms)
        except Exception as exc:  # noqa: BLE001 - degrade to an error column
            entry.update({"error": str(exc), "_structure": None, "_grid": None})
        resolved_targets.append(entry)

    # --- Dock every (target, ligand) pair.
    cells = []
    for target in resolved_targets:
        for ligand in resolved_ligands:
            cell = {"target_id": target["id"], "ligand_id": ligand["id"]}
            if target.get("error"):
                cell["error"] = f"target unavailable: {target['error']}"
            elif ligand.get("error"):
                cell["error"] = f"ligand unavailable: {ligand['error']}"
            else:
                try:
                    pose = docking.dock_with_grid(
                        target["_grid"], ligand["_atoms"],
                        tuple(target["center"]),
                        seeds=seeds, mc_steps=mc_steps, seed=random_seed,
                    )
                    res_name = (ligand.get("formula") or "LIG")[:3].upper()
                    component = docking.pose_to_component(pose, res_name)
                    profile = interactions.profile_component(
                        target["_structure"], component
                    )
                    cell.update({
                        "error": None,
                        "score": pose["score"],
                        "ligand_efficiency": pose["ligand_efficiency"],
                        "n_heavy_atoms": pose["n_heavy_atoms"],
                        **_profile_summary(profile),
                    })
                    if measured_activity is not None:
                        try:
                            cell["measured"] = measured_activity(ligand, target)
                        except Exception:  # noqa: BLE001 - measurement is a bonus
                            cell["measured"] = None
                except Exception as exc:  # noqa: BLE001 - one bad cell only
                    cell["error"] = str(exc)
            cells.append(cell)
            done += 1
            if progress is not None:
                progress(done, total, f"{target['id']} x {ligand['id']}")

    normalization, selectivity = _normalize(cells, resolved_targets, resolved_ligands)

    # Drop the heavy private handles before the result crosses a JSON boundary.
    for entry in resolved_targets:
        entry.pop("_structure", None)
        entry.pop("_grid", None)
    for entry in resolved_ligands:
        entry.pop("_atoms", None)

    scored = [c for c in cells if c.get("error") is None]
    return {
        "targets": resolved_targets,
        "ligands": resolved_ligands,
        "cells": cells,
        "n_cells": len(cells),
        "n_scored": len(scored),
        "n_failed": len(cells) - len(scored),
        "normalization": {
            "basis": "ligand_efficiency (grid score / heavy atom)",
            "per_target": normalization,
            "note": (
                "z_target standardizes ligand efficiency within one target "
                "column (sign flipped: higher = better than this panel's "
                "average for that target). Raw scores are NOT comparable "
                "across targets; only z-scores are. Every z is relative to "
                "the submitted panel and changes if the panel changes."
            ),
        },
        "selectivity": selectivity,
        "methods": methods_block(resolved_targets, resolved_ligands,
                                 seeds, mc_steps, random_seed),
    }


def methods_block(targets, ligands, seeds, mc_steps, random_seed):
    """Reproducibility record for a panel run (mirrors server._methods_block)."""
    return {
        "tool": f"SnaCleX v{__version__}",
        "analysis": "ligand x target interaction panel",
        "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"
        ),
        "n_targets": len(targets),
        "n_ligands": len(ligands),
        "receptor_prep": (
            "Rigid receptor per target; protein heavy atoms only -- waters, "
            "ions, cofactors and co-crystallized ligands are excluded from "
            "every scoring grid. No explicit hydrogens; heavy-atom geometry at "
            "implicit pH ~7. Protonation state is therefore NOT modelled, which "
            "matters for amine ligands and catalytic residues."
        ),
        "box": {
            "edge_A": round(2 * docking.GRID_HALF, 1),
            "grid_spacing_A": docking.SPACING,
            "translation_search_A": docking.TRANS_HALF,
        },
        "search": {
            "algorithm": "Monte-Carlo rigid-body, simulated-annealing acceptance",
            "seeds": seeds,
            "mc_steps": mc_steps,
            "random_seed": random_seed,
            "note": (
                "Panels use a lighter search than a single focused dock "
                f"({PANEL_SEEDS} vs 220 seeds) to keep M x N tractable; "
                "re-dock any pair of interest at full budget before relying "
                "on its pose."
            ),
        },
        "scoring": (
            "AutoDock-style grid-map empirical score: steric (smoothed "
            "Lennard-Jones) + hydrogen-bond + hydrophobic channels, "
            "trilinear-interpolated. Relative units (lower = better) -- NOT "
            "calibrated to kcal/mol and NOT an affinity."
        ),
        "normalization": (
            "Per-target standardization of ligand efficiency before any "
            "cross-target comparison; see result.normalization."
        ),
        "interaction_cutoffs_A": {
            "hydrogen_bond": interactions.HB_MAX,
            "salt_bridge": interactions.SALT_MAX,
            "hydrophobic": interactions.HYDRO_MAX,
            "metal_coordination": interactions.METAL_MAX,
            "aromatic_centroid": interactions.ARO_CENTROID_MAX,
        },
        "limitations": [
            "Rigid ligand (one PubChem conformer) and rigid receptor: induced "
            "fit, side-chain rearrangement and ligand strain are all ignored.",
            "No protonation/tautomer enumeration. Catecholamines and other "
            "amines are charged at physiological pH; the neutral heavy-atom "
            "geometry used here misstates their electrostatics.",
            "Cofactor-dependent enzymes (FAD, SAM, PLP, tetrahydrobiopterin, "
            "metal centers) are stripped from the grid, so substrate sites that "
            "require the cofactor are scored in an unphysical empty state.",
            "Scores are panel-relative contrasts, not affinities, and carry no "
            "concentration, expression level, competition or time dependence.",
            "A panel ranks hypotheses for measurement. It does not predict a "
            "physiological or clinical outcome.",
        ],
        "disclaimer": (
            "Research-only. A panel is a structure-derived hypothesis matrix, "
            "not a simulation of a biological event and not clinical guidance. "
            "Confirm any pair of interest against measured binding data and "
            "orthogonal experiment."
        ),
    }
