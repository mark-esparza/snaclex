"""Worked example for the SnaCleX preprint: trypsin + benzamidine (3PTB)."""
import json, math, time
from snaclex import docking, evolution, interactions, pdbparse, pockets, pubchem, rcsb

out = {"pdb": "3PTB", "run_utc": time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())}
s = pdbparse.parse_structure(rcsb.fetch_structure("3PTB"))
ben = [c for c in s.ligand_components if c.res_name == "BEN"][0]
center = docking.component_center(ben)

prof = interactions.profile_component(s, ben)
out["interactions"] = {"contacts": prof.get("contacts") or prof.get("interactions"),
                       "contact_residues": prof.get("contact_residues"),
                       "keys": list(prof.keys())}

found = pockets.detect_pockets(s)
out["pockets_top5"] = [{"rank": i + 1, "volume": p.get("volume"), "enclosure": p.get("mean_psp") or p.get("enclosure"),
                        "score": p.get("score"), "tier": p.get("tier"),
                        "dist_to_ligand_A": round(math.dist(p["center"], center), 2),
                        "lining": [r.get("res_id") for r in (p.get("lining") or p.get("lining_residues") or [])][:12],
                        "keys": list(p.keys())} for i, p in enumerate(found[:5])]
out["n_pockets"] = len(found)

try:
    evo = evolution.analyze(s, rcsb.fetch_uniprot_accessions("3PTB"))
    if evo and evo.get("available") is not False:
        pc = evolution.annotate_pockets(evo, found)
        out["conservation"] = {"pfam": evo.get("pfam"), "n_seqs": evo.get("n_sequences") or evo.get("n_seqs"),
                               "keys": list(evo.keys()), "pocket_conservation_top3": pc[:3]}
    else:
        out["conservation"] = evo
except Exception as e:
    out["conservation"] = {"error": str(e)}

grid = docking.build_grid(s, center)
cmpd = pubchem.lookup_compound("benzamidine")
lig = pubchem.fetch_3d_atoms(cmpd["cid"])
t0 = time.time()
pose = docking.dock_with_grid(grid, lig["atoms"], center, seed=0)
pc_prof = interactions.profile_component(s, docking.pose_to_component(pose, "BEN"))
out["redock_pubchem"] = {"cid": cmpd["cid"], "score": pose["score"], "rmsd_A": docking.rmsd_to_reference(pose, ben),
                         "seconds": round(time.time() - t0, 1),
                         "contact_residues": pc_prof.get("contact_residues")}

screen = []
for name in ["benzamidine", "4-aminobenzamidine", "benzylamine", "benzoic acid", "phenol", "toluene"]:
    try:
        c = pubchem.lookup_compound(name); a = pubchem.fetch_3d_atoms(c["cid"])["atoms"]
        p = docking.dock_with_grid(grid, a, center, seeds=160, seed=0)
        heavy = sum(1 for x in a if x["element"] != "H")
        screen.append({"name": name, "cid": c["cid"], "heavy": heavy, "score": p["score"], "per_atom": round(p["score"] / heavy, 3)})
    except Exception as e:
        screen.append({"name": name, "error": str(e)})
out["screen"] = sorted(screen, key=lambda r: r.get("score", 1e9))
json.dump(out, open("example_results.json", "w"), indent=1, default=str)
print("DONE example")
