"""Single-process timing on 20 Astex complexes (every 4th id), seed 0, nothing else running."""
import json, os, statistics as st, time
from run_posebusters import ROOT, read_sdf_heavy
from snaclex import benchmark, docking, interactions, pdbparse, pockets

ids = [l.strip() for l in open(os.path.join(ROOT, "astex_diverse_set_ids.txt")) if l.strip()][::4][:20]
rows = []
for cid in ids:
    d = os.path.join(ROOT, "astex_diverse_set", cid)
    t0 = time.perf_counter()
    s = pdbparse.parse_structure(open(os.path.join(d, f"{cid}_protein.pdb"), encoding="utf-8").read())
    comp = read_sdf_heavy(os.path.join(d, f"{cid}_ligand.sdf"), cid.split("_")[1][:3])
    t1 = time.perf_counter(); interactions.profile_component(s, comp)
    t2 = time.perf_counter(); pockets.detect_pockets(s)
    t3 = time.perf_counter(); c = docking.component_center(comp); g = docking.build_grid(s, c)
    t4 = time.perf_counter()
    docking.dock_with_grid(g, [{"element": a.element, "x": a.x, "y": a.y, "z": a.z} for a in comp.atoms], c, seed=0)
    t5 = time.perf_counter()
    rows.append({"id": cid, "atoms": len(s.protein_atoms), "parse": t1 - t0, "interactions": t2 - t1,
                 "pockets": t3 - t2, "grid": t4 - t3, "dock": t5 - t4, "total": t5 - t0})
    print(cid, {k: round(v, 2) for k, v in rows[-1].items() if isinstance(v, float)}, flush=True)
summ = {k: round(st.median(r[k] for r in rows), 2) for k in ("parse", "interactions", "pockets", "grid", "dock", "total")}
summ["median_protein_atoms"] = st.median(r["atoms"] for r in rows)
json.dump({"summary": summ, "rows": rows}, open("timing.json", "w"), indent=1)
print("SUMMARY", summ)
