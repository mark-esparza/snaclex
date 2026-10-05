"""SnaCleX validation on the Astex Diverse (85) and PoseBusters Benchmark (428) sets.

Rigid self-docking of the crystal ligand conformer, 3 random seeds per complex.
Writes one JSON line per complex to <set>_results.jsonl (resumable) and the
docked poses to poses/<set>/<id>_seed<k>.xyz for symmetry-corrected RMSD later.
"""
import json, math, os, sys, time, traceback
from multiprocessing import Pool

from snaclex import benchmark, docking, interactions, pdbparse, pockets

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pb")
SEEDS = (0, 1, 2)


def read_sdf_heavy(path, res_name):
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    n_atoms = int(lines[3][0:3])
    atoms = []
    for i, ln in enumerate(lines[4:4 + n_atoms]):
        x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
        el = ln[31:34].strip().upper()
        if el in ("H", "D"):
            continue
        atoms.append(pdbparse.Atom(serial=i + 1, name=f"{el}{i + 1}", res_name=res_name,
                                   chain="L", res_seq=1, icode="", x=x, y=y, z=z,
                                   element=el, is_hetero=True))
    return pdbparse.Component(res_name, "L", 1, "", atoms)


def matched_rmsd(pose, comp):
    ref = [(a.x, a.y, a.z) for a in comp.atoms]
    s = sum((x - r[0]) ** 2 + (y - r[1]) ** 2 + (z - r[2]) ** 2
            for (x, y, z), r in zip(pose["pose_coords"], ref))
    return round(math.sqrt(s / len(ref)), 3)


def run_case(args):
    set_name, cid = args
    d = os.path.join(ROOT, set_name, cid)
    out = {"set": set_name, "id": cid}
    try:
        t0 = time.time()
        with open(os.path.join(d, f"{cid}_protein.pdb"), encoding="utf-8") as fh:
            s = pdbparse.parse_structure(fh.read())
        comp = read_sdf_heavy(os.path.join(d, f"{cid}_ligand.sdf"), cid.split("_")[1][:3])
        center = docking.component_center(comp)
        out["n_heavy_atoms"] = len(comp.atoms)
        out["n_protein_atoms"] = len(s.protein_atoms)
        tp = time.time()
        found = pockets.detect_pockets(s)
        out["pocket"] = benchmark.pocket_recovery(found, center, 6.0)
        out["pocket_seconds"] = round(time.time() - tp, 2)
        ref_res = {r["res_id"] for r in interactions.profile_component(s, comp)["contact_residues"]}
        out["interactions_total"] = len(ref_res)
        tg = time.time()
        grid = docking.build_grid(s, center)
        out["grid_seconds"] = round(time.time() - tg, 2)
        lig = [{"element": a.element, "x": a.x, "y": a.y, "z": a.z} for a in comp.atoms]
        pose_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "poses", set_name)
        os.makedirs(pose_dir, exist_ok=True)
        runs = []
        for seed in SEEDS:
            td = time.time()
            pose = docking.dock_with_grid(grid, lig, center, seeds=220, mc_steps=40, seed=seed)
            dock_s = time.time() - td
            pc = docking.pose_to_component(pose, comp.res_name)
            pose_res = {r["res_id"] for r in interactions.profile_component(s, pc)["contact_residues"]}
            with open(os.path.join(pose_dir, f"{cid}_seed{seed}.xyz"), "w") as fh:
                fh.write("\n".join(f"{e} {x:.4f} {y:.4f} {z:.4f}" for (x, y, z), e in
                                   zip(pose["pose_coords"], pose["elements"])))
            runs.append({
                "seed": seed,
                "score": pose["score"],
                "rmsd_matched_A": matched_rmsd(pose, comp),
                "rmsd_nearest_A": docking.rmsd_to_reference(pose, comp),
                "contacts_recovered": len(ref_res & pose_res),
                "clash": benchmark._has_clash(pose, s, 2.0),
                "dock_seconds": round(dock_s, 2),
            })
        out["runs"] = runs
        out["total_seconds"] = round(time.time() - t0, 2)
    except Exception:
        out["error"] = traceback.format_exc(limit=2)
    return out


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    only = sys.argv[2] if len(sys.argv) > 2 else None
    for set_name in ("astex_diverse_set", "posebusters_benchmark_set"):
        ids = [l.strip() for l in open(os.path.join(ROOT, f"{set_name}_ids.txt")) if l.strip()]
        if only:
            ids = [i for i in ids if i == only]
            if not ids:
                continue
        res_path = f"{set_name}_results.jsonl"
        done = set()
        if os.path.exists(res_path):
            done = {json.loads(l)["id"] for l in open(res_path) if l.strip()}
        todo = [(set_name, i) for i in ids if i not in done]
        print(f"{set_name}: {len(ids)} total, {len(todo)} to run", flush=True)
        if only:
            print(json.dumps(run_case(todo[0]) if todo else {}, indent=1), flush=True)
            continue
        with Pool(workers) as pool, open(res_path, "a", encoding="utf-8") as fh:
            for k, r in enumerate(pool.imap_unordered(run_case, todo), 1):
                fh.write(json.dumps(r) + "\n"); fh.flush()
                rm = [x["rmsd_matched_A"] for x in r.get("runs", [])]
                print(f"[{k}/{len(todo)}] {r['id']} rmsd={rm} err={'error' in r}", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
