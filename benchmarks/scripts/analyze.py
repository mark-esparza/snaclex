"""Summarize the SnaCleX benchmark and PLIP comparison (run with the plip env python)."""
import json, os, statistics as st, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "pb")
SETS = [s for s in ("astex_diverse_set", "posebusters_benchmark_set")
        if os.path.exists(os.path.join(HERE, f"{s}_results.jsonl"))]


def load(path):
    rows = {}
    for l in open(path, encoding="utf-8"):
        if l.strip():
            r = json.loads(l); rows[r["id"]] = r
    return list(rows.values())


def sym_rmsd(set_name, cid, seed):
    ref = Chem.MolFromMolFile(os.path.join(ROOT, set_name, cid, f"{cid}_ligand.sdf"), sanitize=False)
    if ref is None:
        return None
    ref = Chem.RemoveHs(ref, sanitize=False)
    xyz = [l.split() for l in open(os.path.join(HERE, "poses", set_name, f"{cid}_seed{seed}.xyz"))]
    if len(xyz) != ref.GetNumAtoms():
        return None
    probe = Chem.Mol(ref)
    conf = probe.GetConformer()
    for i, (_, x, y, z) in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    try:
        ref.UpdatePropertyCache(strict=False); probe.UpdatePropertyCache(strict=False)
        return round(rdMolAlign.CalcRMS(probe, ref), 3)
    except Exception:
        return None


def pct(a, b):
    return f"{a}/{b} ({100 * a / b:.1f}%)" if b else "n/a"


def docking_summary(set_name):
    rows = load(os.path.join(HERE, f"{set_name}_results.jsonl"))
    ok = [r for r in rows if "runs" in r]
    errs = [r for r in rows if "error" in r]
    for r in ok:
        for run in r["runs"]:
            run["rmsd_sym_A"] = sym_rmsd(set_name, r["id"], run["seed"])
    n = len(ok)
    out = {"set": set_name, "n_total": len(rows), "n_scored": n, "n_errors": len(errs),
           "errors": [(r["id"], r["error"].strip().splitlines()[-1]) for r in errs]}

    def rm(run):
        return run["rmsd_sym_A"] if run["rmsd_sym_A"] is not None else run["rmsd_matched_A"]

    per_seed = []
    for k in range(3):
        per_seed.append(sum(1 for r in ok if rm(r["runs"][k]) <= 2.0))
    out["success_per_seed"] = per_seed
    out["success_mean_pct"] = round(100 * st.fmean(per_seed) / n, 1)
    out["success_sd_pct"] = round(100 * st.pstdev(per_seed) / n, 1)
    best = [min(r["runs"], key=lambda x: x["score"]) for r in ok]
    out["success_best_score_of_3"] = pct(sum(1 for b in best if rm(b) <= 2.0), n)
    out["success_any_of_3"] = pct(sum(1 for r in ok if min(rm(x) for x in r["runs"]) <= 2.0), n)
    out["median_rmsd_seed0_A"] = round(st.median(rm(r["runs"][0]) for r in ok), 2)
    out["median_rmsd_best_score_A"] = round(st.median(rm(b) for b in best), 2)
    out["nearest_atom_success_seed0"] = pct(sum(1 for r in ok if r["runs"][0]["rmsd_nearest_A"] <= 2.0), n)
    out["matched_vs_sym_seed0"] = pct(sum(1 for r in ok if r["runs"][0]["rmsd_matched_A"] <= 2.0), n)
    out["n_sym_rmsd_fallback"] = sum(1 for r in ok for x in r["runs"] if x["rmsd_sym_A"] is None)
    out["score_rank_correct"] = pct(sum(1 for r in ok if rm(min(r["runs"], key=lambda x: x["score"])) ==
                                        min(rm(x) for x in r["runs"])), n)
    out["pocket_within_6A"] = pct(sum(1 for r in ok if r["pocket"]["recovered"]), n)
    out["pocket_top1"] = pct(sum(1 for r in ok if r["pocket"]["recovered"] and r["pocket"]["rank"] == 1), n)
    out["pocket_top3"] = pct(sum(1 for r in ok if r["pocket"]["recovered"] and r["pocket"]["rank"] <= 3), n)
    out["no_pocket_found"] = sum(1 for r in ok if r["pocket"]["n_pockets"] == 0)
    cr = [x["contacts_recovered"] / r["interactions_total"] for r in ok for x in r["runs"][:1] if r["interactions_total"]]
    out["contact_recovery_mean_seed0"] = round(100 * st.fmean(cr), 1)
    out["clash_free_seed0"] = pct(sum(1 for r in ok if not r["runs"][0]["clash"]), n)
    out["median_seconds_total"] = round(st.median(r["total_seconds"] for r in ok), 1)
    out["median_seconds_dock"] = round(st.median(x["dock_seconds"] for r in ok for x in r["runs"]), 2)
    out["median_seconds_pocket"] = round(st.median(r["pocket_seconds"] for r in ok), 2)
    out["median_heavy_atoms"] = st.median(r["n_heavy_atoms"] for r in ok)
    bins = [(0, 20), (21, 35), (36, 999)]
    out["by_size_seed0"] = {f"{a}-{b}": pct(sum(1 for r in ok if a <= r["n_heavy_atoms"] <= b and rm(r["runs"][0]) <= 2.0),
                                            sum(1 for r in ok if a <= r["n_heavy_atoms"] <= b)) for a, b in bins}
    with open(os.path.join(HERE, f"{set_name}_with_sym.json"), "w") as fh:
        json.dump(ok, fh)
    return out


def plip_summary(set_name):
    path = os.path.join(HERE, f"plip_{set_name}.jsonl")
    if not os.path.exists(path):
        return None
    rows = [r for r in load(path)]
    ok = [r for r in rows if "plip" in r]
    out = {"set": set_name, "n_total": len(rows), "n_compared": len(ok),
           "errors": [(r["id"], r["error"].strip().splitlines()[-1]) for r in rows if "error" in r][:10]}
    types = ("hydrogen_bond", "salt_bridge", "hydrophobic", "aromatic", "metal")
    for t in types + ("any",):
        tp = fp = fn = 0
        for r in ok:
            if t == "any":
                p = set().union(*[set(r["plip"][k]) for k in types + ("halogen",)])
                s = set().union(*[set(r["snaclex"][k]) for k in types])
            else:
                p, s = set(r["plip"][t]), set(r["snaclex"][t])
            tp += len(p & s); fp += len(s - p); fn += len(p - s)
        prec = tp / (tp + fp) if tp + fp else None
        rec = tp / (tp + fn) if tp + fn else None
        f1 = 2 * prec * rec / (prec + rec) if prec and rec else None
        out[t] = {"plip_residues": tp + fn, "snaclex_residues": tp + fp, "shared": tp,
                  "precision": round(prec, 3) if prec is not None else None,
                  "recall": round(rec, 3) if rec is not None else None,
                  "f1": round(f1, 3) if f1 is not None else None}
    jac = []
    for r in ok:
        p = set().union(*[set(r["plip"][k]) for k in types + ("halogen",)])
        s = set().union(*[set(r["snaclex"][k]) for k in types])
        if p | s:
            jac.append(len(p & s) / len(p | s))
    out["any_jaccard_median"] = round(st.median(jac), 3) if jac else None
    out["any_jaccard_mean"] = round(st.fmean(jac), 3) if jac else None
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = {}
    for s in SETS:
        if which in ("all", "dock"):
            res[f"dock_{s}"] = docking_summary(s)
    for s in ("astex_diverse_set", "posebusters_benchmark_set"):
        if which in ("all", "plip"):
            p = plip_summary(s)
            if p:
                res[f"plip_{s}"] = p
    json.dump(res, open(os.path.join(HERE, "summary.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))
