"""Compute SnaCleX outputs for 3PTB and write three 3Dmol.js figure pages."""
import json, math
from snaclex import docking, evolution, interactions, pdbparse, pockets, pubchem, rcsb

pdb_text = rcsb.fetch_structure("3PTB")
s = pdbparse.parse_structure(pdb_text)
ben = [c for c in s.ligand_components if c.res_name == "BEN"][0]
center = docking.component_center(ben)
prof = interactions.profile_component(s, ben)
found = pockets.detect_pockets(s)
cons = {}
try:
    evo = evolution.analyze(s, rcsb.fetch_uniprot_accessions("3PTB"))
    for a in evolution.annotate_pockets(evo, found):
        cons[a["index"]] = a["mean_conservation"]
except Exception as e:
    print("conservation failed", e)
grid = docking.build_grid(s, center)
cid = pubchem.lookup_compound("benzamidine")["cid"]
lig = pubchem.fetch_3d_atoms(cid)["atoms"]
pose = docking.dock_with_grid(grid, lig, center, seed=0)
rmsd = docking.rmsd_to_reference(pose, ben)
print("redock rmsd", rmsd, "score", pose["score"])

data = {
    "pdb": pdb_text,
    "pose_pdb": docking.pose_to_pdb(pose, "DCK"),
    "rmsd": rmsd,
    "contacts": [{"type": r["type"], "l": r["ligand_atom"]["xyz"], "p": r["protein_atom"]["xyz"]}
                 for r in prof["interactions"]],
    "residues": [{"resi": r["res_seq"], "resn": r["res_name"], "chain": r["chain"]}
                 for r in prof["contact_residues"]],
    "pockets": [{"rank": p["index"] + 1, "center": p["center"], "points": p["points"],
                 "volume": p["volume_A3"], "score": p["score"], "cons": cons.get(p["index"]),
                 "dist": round(math.dist(p["center"], center), 1)} for p in found],
}
json.dump({k: v for k, v in data.items() if k not in ("pdb",)},
          open("fig_data.json", "w"), indent=1)
tmpl = open("fig_template.html", encoding="utf-8").read()
for fig in ("fig1", "fig2", "fig3"):
    html = tmpl.replace("__DATA__", json.dumps(data)).replace("__FIG__", fig)
    open(f"{fig}.html", "w", encoding="utf-8").write(html)
print("pages written")
