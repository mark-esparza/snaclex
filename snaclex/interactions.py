"""Atomic-level protein-ligand interaction profiler.

Given a parsed Structure and a chosen hetero Component (ligand, ion, or metal),
this computes heavy-atom contacts between the component and the protein and
classifies each into an interaction type using geometric criteria. This mirrors
the approach of tools like PLIP, simplified to run dependency-free on
experimental coordinates (no explicit hydrogens assumed).

Criteria are research heuristics, not force-field calculations.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .pdbparse import Component, Structure

# Distance cutoffs in Angstroms.
HB_MAX = 3.6          # N/O .. N/O polar contact (heavy-atom)
HB_S_MAX = 3.9        # contacts involving sulfur
SALT_MAX = 4.0        # charged group .. charged group
HYDRO_MIN = 2.8       # avoid covalent-bond-like distances
HYDRO_MAX = 4.0       # C .. C hydrophobic
METAL_MAX = 2.9       # metal .. coordinating atom
ARO_CENTROID_MAX = 4.5  # aromatic ring centroid .. nearest ligand heavy atom

POLAR_ELEMENTS = {"N", "O"}

# Charged protein side-chain atoms.
PROT_POS_ATOMS = {
    ("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"),
    ("LYS", "NZ"),
    ("HIS", "ND1"), ("HIS", "NE2"),
}
PROT_NEG_ATOMS = {
    ("ASP", "OD1"), ("ASP", "OD2"),
    ("GLU", "OE1"), ("GLU", "OE2"),
}

# Aromatic ring atoms by residue, for centroid-based pi detection.
AROMATIC_RINGS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

INTERACTION_TYPES = [
    "metal_coordination",
    "salt_bridge",
    "hydrogen_bond",
    "hydrophobic",
    "aromatic",
]


@dataclass
class _Grid:
    cell: float
    buckets: dict

    @classmethod
    def build(cls, atoms, cell: float):
        buckets: dict = {}
        for a in atoms:
            key = (int(a.x // cell), int(a.y // cell), int(a.z // cell))
            buckets.setdefault(key, []).append(a)
        return cls(cell, buckets)

    def neighbors(self, x: float, y: float, z: float):
        cx, cy, cz = int(x // self.cell), int(y // self.cell), int(z // self.cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    bucket = self.buckets.get((cx + dx, cy + dy, cz + dz))
                    if bucket:
                        yield from bucket


def _dist(a, b) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def _atom_json(a, include_res: bool):
    d = {
        "name": a.name,
        "element": a.element,
        "serial": a.serial,
        "xyz": [round(a.x, 3), round(a.y, 3), round(a.z, 3)],
    }
    if include_res:
        d.update(
            res_name=a.res_name,
            res_seq=a.res_seq,
            chain=a.chain,
            res_id=f"{a.chain}/{a.res_name}{a.res_seq}",
        )
    return d


def _classify(lig, prot, d: float, is_metal_comp: bool) -> str | None:
    le, pe = lig.element, prot.element
    pkey = (prot.res_name, prot.name)

    # Metal coordination (component is a metal ion).
    if is_metal_comp and le in {"NA", "K", "MG", "CA", "MN", "FE", "CO", "NI",
                                "CU", "ZN", "MO", "W", "CD", "HG", "PT", "AU",
                                "AG", "PB", "BA", "SR", "LI", "AL", "V", "CR"}:
        if pe in {"O", "N", "S"} and d <= METAL_MAX:
            return "metal_coordination"
        return None

    # Salt bridge (charged .. charged), heuristic on protonation.
    if d <= SALT_MAX:
        if pkey in PROT_NEG_ATOMS and le == "N":
            return "salt_bridge"
        if pkey in PROT_POS_ATOMS and le == "O":
            return "salt_bridge"

    # Hydrogen bond / polar contact (heavy atom).
    if le in POLAR_ELEMENTS and pe in POLAR_ELEMENTS and d <= HB_MAX:
        return "hydrogen_bond"
    if ("S" in (le, pe)) and (le in POLAR_ELEMENTS or pe in POLAR_ELEMENTS) and d <= HB_S_MAX:
        return "hydrogen_bond"

    # Hydrophobic carbon contact.
    if le == "C" and pe == "C" and HYDRO_MIN <= d <= HYDRO_MAX:
        return "hydrophobic"

    return None


def _aromatic_interactions(component: Component, structure: Structure):
    """Detect aromatic ring centroid .. ligand contacts (possible pi-stacking)."""
    lig_heavy = [a for a in component.atoms if a.element != "H"]
    if not lig_heavy:
        return []

    # Group protein atoms by residue to assemble ring atom sets.
    by_res: dict[tuple, list] = {}
    for a in structure.protein_atoms:
        if a.res_name in AROMATIC_RINGS:
            by_res.setdefault((a.chain, a.res_seq, a.res_name), []).append(a)

    results = []
    for (chain, res_seq, res_name), res_atoms in by_res.items():
        ring_names = set(AROMATIC_RINGS[res_name])
        ring_atoms = [a for a in res_atoms if a.name in ring_names]
        if len(ring_atoms) < 3:
            continue
        cx = sum(a.x for a in ring_atoms) / len(ring_atoms)
        cy = sum(a.y for a in ring_atoms) / len(ring_atoms)
        cz = sum(a.z for a in ring_atoms) / len(ring_atoms)

        nearest = None
        nearest_d = ARO_CENTROID_MAX
        for la in lig_heavy:
            dd = math.sqrt((la.x - cx) ** 2 + (la.y - cy) ** 2 + (la.z - cz) ** 2)
            if dd < nearest_d:
                nearest_d = dd
                nearest = la
        if nearest is not None:
            centroid_atom = type(ring_atoms[0])(
                serial=ring_atoms[0].serial, name="ring-centroid",
                res_name=res_name, chain=chain, res_seq=res_seq, icode="",
                x=cx, y=cy, z=cz, element="C", is_hetero=False,
            )
            results.append(
                {
                    "type": "aromatic",
                    "distance": round(nearest_d, 2),
                    "ligand_atom": _atom_json(nearest, include_res=False),
                    "protein_atom": _atom_json(centroid_atom, include_res=True),
                }
            )
    return results


# Protein-protein interface (paratope/epitope) analysis --------------------

INTERFACE_MAX = 4.5   # heavy-atom contact distance across a protein interface
_PROBE = 1.4          # water-probe radius (Angstrom) for the SASA estimate
# Coarse per-element van der Waals radii (heavy atoms only; no hydrogens).
_VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
# A small Fibonacci sphere for Shrake-Rupley; 96 points balances speed/accuracy.
_SASA_POINTS = 96


def _sphere_points(n: int):
    import math as _m
    pts = []
    phi = _m.pi * (3.0 - _m.sqrt(5.0))
    for i in range(n):
        y = 1.0 - 2.0 * (i + 0.5) / n
        r = _m.sqrt(max(0.0, 1.0 - y * y))
        theta = phi * i
        pts.append((_m.cos(theta) * r, y, _m.sin(theta) * r))
    return pts


def _residue_sasa(atoms, neighbor_atoms):
    """Shrake-Rupley SASA per residue over ``atoms`` (Å²), occluded by neighbours.

    ``neighbor_atoms`` supplies the occluding context (may include atoms outside
    ``atoms``, e.g. the partner chain when computing bound-state SASA).
    """
    heavy = [a for a in atoms if a.element in _VDW]
    if not heavy:
        return {}
    occl = [a for a in neighbor_atoms if a.element in _VDW]
    grid = _Grid.build(occl, cell=8.0)
    points = _sphere_points(_SASA_POINTS)
    per_res: dict[tuple, float] = {}
    for a in heavy:
        ra = _VDW[a.element] + _PROBE
        near = [b for b in grid.neighbors(a.x, a.y, a.z)
                if b is not a and _dist(a, b) < ra + _VDW.get(b.element, 1.7) + _PROBE]
        accessible = 0
        for px, py, pz in points:
            sx, sy, sz = a.x + px * ra, a.y + py * ra, a.z + pz * ra
            buried = False
            for b in near:
                rb = _VDW.get(b.element, 1.7) + _PROBE
                if (sx - b.x) ** 2 + (sy - b.y) ** 2 + (sz - b.z) ** 2 < rb * rb:
                    buried = True
                    break
            if not buried:
                accessible += 1
        area = 4.0 * math.pi * ra * ra * accessible / _SASA_POINTS
        per_res[(a.chain, a.res_seq)] = per_res.get((a.chain, a.res_seq), 0.0) + area
    return per_res


def profile_interface(structure: Structure, ab_chains, ag_chains) -> dict:
    """Analyze a protein-protein (antibody paratope / antigen epitope) interface.

    Reuses the atomic contact grid — this is a conditional branch of the
    interaction logic, not a new engine. Returns the paratope residues (on the
    antibody side), the epitope residues (on the antigen side), and a coarse
    buried-surface-area estimate. BSA uses a Shrake-Rupley SASA (heavy atoms, no
    hydrogens) restricted to interface-proximal residues, so it is an estimate.
    """
    ab_set, ag_set = set(ab_chains), set(ag_chains)
    ab_atoms = [a for a in structure.protein_atoms if a.chain in ab_set]
    ag_atoms = [a for a in structure.protein_atoms if a.chain in ag_set]
    if not ab_atoms or not ag_atoms:
        return {"available": False, "reason": "Missing antibody or antigen chain atoms."}

    ag_grid = _Grid.build(ag_atoms, cell=max(6.0, INTERFACE_MAX + 1.0))
    paratope: dict[tuple, dict] = {}
    epitope: dict[tuple, dict] = {}
    n_contacts = 0
    for a in ab_atoms:
        for b in ag_grid.neighbors(a.x, a.y, a.z):
            d = _dist(a, b)
            if d > INTERFACE_MAX:
                continue
            n_contacts += 1
            for atom, store in ((a, paratope), (b, epitope)):
                key = (atom.chain, atom.res_seq)
                e = store.setdefault(key, {
                    "chain": atom.chain, "res_name": atom.res_name,
                    "res_seq": atom.res_seq, "contacts": 0, "min_distance": d,
                })
                e["contacts"] += 1
                e["min_distance"] = min(e["min_distance"], round(d, 2))

    def _residue_list(store):
        out = sorted(store.values(), key=lambda e: (-e["contacts"], e["min_distance"]))
        for e in out:
            e["min_distance"] = round(e["min_distance"], 2)
        return out

    # Buried surface area over interface-proximal residues only (bounded).
    iface_keys = set(paratope) | set(epitope)
    bsa = None
    if iface_keys:
        iface_atoms = [a for a in structure.protein_atoms
                       if (a.chain, a.res_seq) in iface_keys]
        ab_iface = [a for a in iface_atoms if a.chain in ab_set]
        ag_iface = [a for a in iface_atoms if a.chain in ag_set]
        # Free-state SASA (each side alone) minus bound-state SASA (occluded by
        # the partner); summed over both sides = total buried area.
        free_ab = _residue_sasa(ab_iface, ab_atoms)
        free_ag = _residue_sasa(ag_iface, ag_atoms)
        bound_ab = _residue_sasa(ab_iface, ab_atoms + ag_atoms)
        bound_ag = _residue_sasa(ag_iface, ag_atoms + ab_atoms)
        buried = sum(max(0.0, free_ab[k] - bound_ab.get(k, 0.0)) for k in free_ab)
        buried += sum(max(0.0, free_ag[k] - bound_ag.get(k, 0.0)) for k in free_ag)
        bsa = round(buried, 1)

    return {
        "available": True,
        "antibody_chains": sorted(ab_set),
        "antigen_chains": sorted(ag_set),
        "atom_contacts": n_contacts,
        "distance_cutoff_A": INTERFACE_MAX,
        "paratope": _residue_list(paratope),
        "epitope": _residue_list(epitope),
        "paratope_residue_count": len(paratope),
        "epitope_residue_count": len(epitope),
        "buried_surface_area_A2": bsa,
    }


def profile_component(structure: Structure, component: Component) -> dict:
    """Return the full atomic interaction profile for one hetero component."""
    is_metal_comp = component.kind in ("metal", "ion")
    grid = _Grid.build(structure.protein_atoms, cell=5.5)

    raw: list[dict] = []
    # For hydrophobic noise control we keep only the closest per (residue).
    hydro_best: dict[tuple, dict] = {}

    for lig in component.atoms:
        if lig.element == "H":
            continue
        for prot in grid.neighbors(lig.x, lig.y, lig.z):
            d = _dist(lig, prot)
            kind = _classify(lig, prot, d, is_metal_comp)
            if kind is None:
                continue
            record = {
                "type": kind,
                "distance": round(d, 2),
                "ligand_atom": _atom_json(lig, include_res=False),
                "protein_atom": _atom_json(prot, include_res=True),
            }
            if kind == "hydrophobic":
                rkey = (prot.chain, prot.res_seq)
                cur = hydro_best.get(rkey)
                if cur is None or d < cur["distance"]:
                    hydro_best[rkey] = record
            else:
                raw.append(record)

    raw.extend(hydro_best.values())
    raw.extend(_aromatic_interactions(component, structure))

    # Counts + per-residue summary.
    counts = {t: 0 for t in INTERACTION_TYPES}
    res_summary: dict[str, dict] = {}
    for r in raw:
        counts[r["type"]] += 1
        pa = r["protein_atom"]
        rid = pa["res_id"]
        entry = res_summary.setdefault(
            rid,
            {
                "res_id": rid,
                "res_name": pa["res_name"],
                "res_seq": pa["res_seq"],
                "chain": pa["chain"],
                "types": set(),
                "total": 0,
                "min_distance": r["distance"],
            },
        )
        entry["types"].add(r["type"])
        entry["total"] += 1
        entry["min_distance"] = min(entry["min_distance"], r["distance"])

    residues = []
    for e in res_summary.values():
        e["types"] = sorted(e["types"])
        residues.append(e)
    residues.sort(key=lambda e: (-e["total"], e["min_distance"]))

    raw.sort(key=lambda r: (r["type"], r["distance"]))

    return {
        "component": {
            "label": component.label,
            "res_name": component.res_name,
            "chain": component.chain,
            "res_seq": component.res_seq,
            "kind": component.kind,
            "atom_count": len(component.atoms),
        },
        "interactions": raw,
        "counts": counts,
        "interaction_total": len(raw),
        "contact_residues": residues,
        "contact_residue_count": len(residues),
    }
