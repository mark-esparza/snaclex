"""Compare SnaCleX interaction profiles with PLIP on crystal complexes.

Run with the 'plip' conda env python. For each complex the PoseBusters protein
PDB and crystal ligand SDF are merged into one PDB (ligand as LIG Z999), PLIP
analyses it, and SnaCleX profiles the same protein + ligand. Residue-level
interaction sets are compared per type. Output: plip_results.jsonl (resumable).
"""
import json, os, sys, tempfile, traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "pb")
sys.path.insert(0, HERE)
from run_posebusters import read_sdf_heavy  # noqa: E402

TYPES = ("hydrogen_bond", "salt_bridge", "hydrophobic", "aromatic", "metal")


def complex_pdb(protein_path, sdf_path):
    prot = [l for l in open(protein_path, encoding="utf-8").read().splitlines()
            if l.startswith(("ATOM", "HETATM"))]
    last = 0
    for l in prot:
        try:
            last = max(last, int(l[6:11]))
        except ValueError:
            pass
    lines = open(sdf_path, encoding="utf-8").read().splitlines()
    na, nb = int(lines[3][0:3]), int(lines[3][3:6])
    out = list(prot) + ["TER"]
    el_n = {}
    for i, ln in enumerate(lines[4:4 + na]):
        x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
        el = ln[31:34].strip()
        el_n[el] = el_n.get(el, 0) + 1
        name = f"{el.upper()}{el_n[el]}"[:4]
        out.append(f"HETATM{last + 1 + i:>5} {name:<4} LIG Z 999    "
                   f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el.upper():>2}")
    con = {}
    for ln in lines[4 + na:4 + na + nb]:
        a, b, order = int(ln[0:3]), int(ln[3:6]), int(ln[6:9])
        for _ in range(order if order in (1, 2, 3) else 1):
            con.setdefault(a, []).append(b); con.setdefault(b, []).append(a)
    for a in sorted(con):
        nbrs = con[a]
        for k in range(0, len(nbrs), 4):
            out.append("CONECT" + f"{last + a:>5}" + "".join(f"{last + b:>5}" for b in nbrs[k:k + 4]))
    out.append("END")
    return "\n".join(out) + "\n"


def plip_sets(pdb_path):
    from plip.structure.preparation import PDBComplex
    mol = PDBComplex()
    mol.output_path = tempfile.gettempdir()
    mol.load_pdb(pdb_path)
    mol.analyze()
    site = next((v for k, v in mol.interaction_sets.items() if k.startswith("LIG:Z:999")), None)
    if site is None:
        raise RuntimeError("PLIP found no LIG site; sites=" + ",".join(mol.interaction_sets))
    key = lambda i: f"{i.reschain}/{i.resnr}"
    s = {t: set() for t in TYPES}
    s["halogen"] = set(); s["water_bridge"] = set()
    s["hydrogen_bond"] = {key(i) for i in site.hbonds_pdon + site.hbonds_ldon}
    s["salt_bridge"] = {key(i) for i in site.saltbridge_lneg + site.saltbridge_pneg}
    s["hydrophobic"] = {key(i) for i in site.hydrophobic_contacts}
    s["aromatic"] = {key(i) for i in site.pistacking + site.pication_laro + site.pication_paro}
    s["metal"] = {key(i) for i in site.metal_complexes}
    s["halogen"] = {key(i) for i in site.halogen_bonds}
    s["water_bridge"] = {key(i) for i in site.water_bridges}
    return s


def snaclex_sets(protein_path, sdf_path, cid):
    from snaclex import interactions, pdbparse
    st = pdbparse.parse_structure(open(protein_path, encoding="utf-8").read())
    comp = read_sdf_heavy(sdf_path, cid.split("_")[1][:3])
    prof = interactions.profile_component(st, comp)
    s = {t: set() for t in TYPES}
    for r in prof["interactions"]:
        pa = r["protein_atom"]
        t = "metal" if r["type"] == "metal_coordination" else r["type"]
        s[t].add(f"{pa['chain']}/{pa['res_seq']}")
    return s


def run_case(args):
    set_name, cid = args
    d = os.path.join(ROOT, set_name, cid)
    prot, sdf = os.path.join(d, f"{cid}_protein.pdb"), os.path.join(d, f"{cid}_ligand.sdf")
    out = {"set": set_name, "id": cid}
    try:
        tmp = os.path.join(tempfile.gettempdir(), f"plipcx_{cid}.pdb")
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write(complex_pdb(prot, sdf))
        p = plip_sets(tmp)
        s = snaclex_sets(prot, sdf, cid)
        os.remove(tmp)
        out["plip"] = {k: sorted(v) for k, v in p.items()}
        out["snaclex"] = {k: sorted(v) for k, v in s.items()}
    except Exception:
        out["error"] = traceback.format_exc(limit=3)
    return out


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    only = sys.argv[2] if len(sys.argv) > 2 else None
    for set_name in ("astex_diverse_set", "posebusters_benchmark_set"):
        ids = [l.strip() for l in open(os.path.join(ROOT, f"{set_name}_ids.txt")) if l.strip()]
        if only:
            if only in ids:
                print(json.dumps(run_case((set_name, only)), indent=1))
            continue
        path = f"plip_{set_name}.jsonl"
        done = {json.loads(l)["id"] for l in open(path)} if os.path.exists(path) else set()
        todo = [(set_name, i) for i in ids if i not in done]
        print(f"{set_name}: {len(todo)} to run", flush=True)
        with Pool(workers) as pool, open(path, "a", encoding="utf-8") as fh:
            for k, r in enumerate(pool.imap_unordered(run_case, todo), 1):
                fh.write(json.dumps(r) + "\n"); fh.flush()
                print(f"[{k}/{len(todo)}] {r['id']} err={'error' in r}", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
