"""Sequential, single-process timing of each SnaCleX stage on the first 20 Astex complexes."""
import json, os, platform, statistics as st, time
from snaclex import docking, interactions, pdbparse, pockets
from run_posebusters import ROOT, read_sdf_heavy

ids = [l.strip() for l in open(os.path.join(ROOT, "astex_diverse_set_ids.txt")) if l.strip()][:20]
rows = []
for cid in ids:
    d = os.path.join(ROOT, "astex_diverse_set", cid)
    t = time.perf_counter()
    s = pdbparse.parse_structure(open(os.path.join(d, f"{cid}_protein.pdb"), encoding="utf-8").read())
    comp = read_sdf_heavy(os.path.join(d, f"{cid}_ligand.sdf"), cid.split("_")[1][:3])
    r = {"id": cid, "protein_atoms": len(s.protein_atoms), "parse": time.perf_counter() - t}
    t = time.perf_counter(); interactions.profile_component(s, comp); r["interactions"] = time.perf_counter() - t
    t = time.perf_counter(); pockets.detect_pockets(s); r["pockets"] = time.perf_counter() - t
    center = docking.component_center(comp)
    t = time.perf_counter(); grid = docking.build_grid(s, center); r["grid"] = time.perf_counter() - t
    lig = [{"element": a.element, "x": a.x, "y": a.y, "z": a.z} for a in comp.atoms]
    t = time.perf_counter(); docking.dock_with_grid(grid, lig, center, seeds=220, mc_steps=40, seed=0); r["dock"] = time.perf_counter() - t
    r["total"] = sum(r[k] for k in ("parse", "interactions", "pockets", "grid", "dock"))
    rows.append(r); print(cid, {k: round(v, 2) for k, v in r.items() if isinstance(v, float)}, flush=True)
summ = {k: round(st.median(r[k] for r in rows), 2) for k in ("parse", "interactions", "pockets", "grid", "dock", "total")}
summ["median_protein_atoms"] = st.median(r["protein_atoms"] for r in rows)
out = {"machine": platform.processor(), "python": platform.python_version(), "n": len(rows), "summary": summ, "rows": rows}
json.dump(out, open("timing.json", "w"), indent=1)
print("SUMMARY", summ)
