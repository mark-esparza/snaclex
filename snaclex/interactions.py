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

from .pdbparse import Atom, Component, Structure

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


# ---------------------------------------------------------------------------
# Pairwise interface profiling (protein-protein / peptide and protein-nucleic).
#
# Both reuse the same spatial grid and the same per-residue summary as the
# ligand profiler above; only the per-contact classification differs. This is
# the substrate for antibody-antigen epitope/paratope mapping and HLA
# peptide-groove analysis (the peptide is just a short chain).
# ---------------------------------------------------------------------------

# Interface contact types (no metal coordination between two polymers).
INTERFACE_TYPES = ["salt_bridge", "hydrogen_bond", "hydrophobic", "aromatic"]

# Nucleotide phosphate-backbone oxygens (carry the negative charge). Includes
# both modern (prime) and legacy (star) atom-name spellings.
NUCLEIC_PHOSPHATE_ATOMS = {
    "P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P",
    "O5'", "O3'", "O5*", "O3*",
}

# Base-ring atoms for stacking detection (centroid proxy). Purines fuse two
# rings; using all ring atoms places the centroid near the fusion, which is an
# adequate proxy for a stacking contact.
_PURINE_RING = ["N9", "C8", "N7", "C5", "C4", "N3", "C2", "N1", "C6"]
_PYRIMIDINE_RING = ["N1", "C2", "N3", "C4", "C5", "C6"]
NUCLEIC_BASE_RINGS = {
    "DA": _PURINE_RING, "DG": _PURINE_RING, "DI": _PURINE_RING,
    "A": _PURINE_RING, "G": _PURINE_RING, "I": _PURINE_RING,
    "DC": _PYRIMIDINE_RING, "DT": _PYRIMIDINE_RING, "DU": _PYRIMIDINE_RING,
    "C": _PYRIMIDINE_RING, "T": _PYRIMIDINE_RING, "U": _PYRIMIDINE_RING,
}

# Ring-centroid .. ring-centroid stacking distance (pi-stacking / T-shaped).
STACK_CENTROID_MAX = 6.0


def _classify_protein_protein(a, b, d: float) -> str | None:
    """Classify a heavy-atom contact between two protein chains."""
    ea, eb = a.element, b.element
    ka, kb = (a.res_name, a.name), (b.res_name, b.name)

    if d <= SALT_MAX:
        if (ka in PROT_POS_ATOMS and kb in PROT_NEG_ATOMS) or (
            ka in PROT_NEG_ATOMS and kb in PROT_POS_ATOMS
        ):
            return "salt_bridge"

    if ea in POLAR_ELEMENTS and eb in POLAR_ELEMENTS and d <= HB_MAX:
        return "hydrogen_bond"
    if ("S" in (ea, eb)) and (ea in POLAR_ELEMENTS or eb in POLAR_ELEMENTS) and d <= HB_S_MAX:
        return "hydrogen_bond"

    if ea == "C" and eb == "C" and HYDRO_MIN <= d <= HYDRO_MAX:
        return "hydrophobic"
    return None


def _classify_protein_nucleic(prot, nuc, d: float) -> str | None:
    """Classify a heavy-atom contact between a protein atom and a nucleic atom."""
    pe, ne = prot.element, nuc.element
    pkey = (prot.res_name, prot.name)

    # Salt bridge: basic side chain .. phosphate backbone.
    if d <= SALT_MAX and pkey in PROT_POS_ATOMS and nuc.name in NUCLEIC_PHOSPHATE_ATOMS:
        return "salt_bridge"

    if pe in POLAR_ELEMENTS and ne in POLAR_ELEMENTS and d <= HB_MAX:
        return "hydrogen_bond"
    if ("S" in (pe, ne)) and (pe in POLAR_ELEMENTS or ne in POLAR_ELEMENTS) and d <= HB_S_MAX:
        return "hydrogen_bond"

    if pe == "C" and ne == "C" and HYDRO_MIN <= d <= HYDRO_MAX:
        return "hydrophobic"
    return None


def _ring_centroids(atoms, ring_map) -> list:
    """Build one ring-centroid pseudo-atom per residue that has a full ring."""
    by_res: dict[tuple, list] = {}
    for a in atoms:
        if a.res_name in ring_map:
            by_res.setdefault((a.chain, a.res_seq, a.res_name), []).append(a)
    centroids = []
    for (chain, res_seq, res_name), res_atoms in by_res.items():
        ring_names = set(ring_map[res_name])
        ring_atoms = [a for a in res_atoms if a.name in ring_names]
        if len(ring_atoms) < 3:
            continue
        n = len(ring_atoms)
        centroids.append(
            Atom(
                serial=ring_atoms[0].serial, name="ring-centroid",
                res_name=res_name, chain=chain, res_seq=res_seq, icode="",
                x=sum(a.x for a in ring_atoms) / n,
                y=sum(a.y for a in ring_atoms) / n,
                z=sum(a.z for a in ring_atoms) / n,
                element="C", is_hetero=False,
            )
        )
    return centroids


def _ring_stacking(atoms_a, rings_a, atoms_b, rings_b, cutoff: float) -> list:
    """Ring-centroid .. ring-centroid contacts (nearest B-ring per A-ring)."""
    ca = _ring_centroids(atoms_a, rings_a)
    cb = _ring_centroids(atoms_b, rings_b)
    out = []
    for x in ca:
        best, best_d = None, cutoff
        for y in cb:
            dd = _dist(x, y)
            if dd < best_d:
                best_d, best = dd, y
        if best is not None:
            out.append(
                {
                    "type": "aromatic",
                    "distance": round(best_d, 2),
                    "atom_a": _atom_json(x, include_res=True),
                    "atom_b": _atom_json(best, include_res=True),
                }
            )
    return out


def _accumulate_residue(store: dict, atom_dict: dict, kind: str, dist: float):
    rid = atom_dict["res_id"]
    entry = store.setdefault(
        rid,
        {
            "res_id": rid, "res_name": atom_dict["res_name"],
            "res_seq": atom_dict["res_seq"], "chain": atom_dict["chain"],
            "types": set(), "total": 0, "min_distance": dist,
        },
    )
    entry["types"].add(kind)
    entry["total"] += 1
    entry["min_distance"] = min(entry["min_distance"], dist)


def _finalize_residues(store: dict) -> list:
    out = []
    for e in store.values():
        e["types"] = sorted(e["types"])
        out.append(e)
    out.sort(key=lambda e: (-e["total"], e["min_distance"]))
    return out


def _profile_pairwise(atoms_a, atoms_b, classify, extra=None) -> dict:
    """Core grid-based contact profile between two atom groups.

    ``classify(a, b, d)`` receives an atom from group A, an atom from group B,
    and their distance, and returns a contact type or None. ``extra`` is a list
    of pre-computed contacts (e.g. ring stacking) to merge in.
    """
    grid = _Grid.build(atoms_b, cell=5.5)
    raw: list[dict] = []
    hydro_best: dict[tuple, dict] = {}

    for a in atoms_a:
        if a.element == "H":
            continue
        for b in grid.neighbors(a.x, a.y, a.z):
            if b.element == "H":
                continue
            d = _dist(a, b)
            kind = classify(a, b, d)
            if kind is None:
                continue
            record = {
                "type": kind,
                "distance": round(d, 2),
                "atom_a": _atom_json(a, include_res=True),
                "atom_b": _atom_json(b, include_res=True),
            }
            if kind == "hydrophobic":
                rkey = (a.chain, a.res_seq, b.chain, b.res_seq)
                cur = hydro_best.get(rkey)
                if cur is None or d < cur["distance"]:
                    hydro_best[rkey] = record
            else:
                raw.append(record)

    raw.extend(hydro_best.values())
    if extra:
        raw.extend(extra)

    counts = {t: 0 for t in INTERFACE_TYPES}
    res_a: dict[str, dict] = {}
    res_b: dict[str, dict] = {}
    for r in raw:
        counts[r["type"]] = counts.get(r["type"], 0) + 1
        _accumulate_residue(res_a, r["atom_a"], r["type"], r["distance"])
        _accumulate_residue(res_b, r["atom_b"], r["type"], r["distance"])

    raw.sort(key=lambda r: (r["type"], r["distance"]))
    return {
        "interactions": raw,
        "counts": counts,
        "interaction_total": len(raw),
        "residues_a": _finalize_residues(res_a),
        "residues_b": _finalize_residues(res_b),
    }


def profile_interface(structure: Structure, chains_a, chains_b) -> dict:
    """Profile the interface between two groups of protein chains.

    ``residues_a`` are the side-A contact residues (e.g. antibody paratope or
    HLA groove); ``residues_b`` are side-B (antigen epitope or bound peptide).
    """
    set_a, set_b = set(chains_a), set(chains_b)
    atoms_a = [a for a in structure.protein_atoms if a.chain in set_a]
    atoms_b = [a for a in structure.protein_atoms if a.chain in set_b]
    if not atoms_a or not atoms_b:
        core = _profile_pairwise([], [], _classify_protein_protein)
    else:
        stacking = _ring_stacking(
            atoms_a, AROMATIC_RINGS, atoms_b, AROMATIC_RINGS, STACK_CENTROID_MAX
        )
        core = _profile_pairwise(
            atoms_a, atoms_b, _classify_protein_protein, extra=stacking
        )
    return {
        "mode": "protein-protein",
        "chains_a": sorted(set_a),
        "chains_b": sorted(set_b),
        "interactions": core["interactions"],
        "counts": core["counts"],
        "interaction_total": core["interaction_total"],
        "interface_residues_a": core["residues_a"],
        "interface_residues_b": core["residues_b"],
        "interface_residue_count": len(core["residues_a"]) + len(core["residues_b"]),
    }


def profile_nucleic_interface(structure: Structure, chain: str | None = None) -> dict:
    """Profile contacts between protein chains and nucleic-acid chains.

    ``protein_residues`` are the DNA/RNA-binding residues; ``nucleic_residues``
    are the contacted nucleotides. Pass ``chain`` to restrict to one nucleic chain.
    """
    prot = structure.protein_atoms
    nuc = [a for a in structure.nucleic_atoms if chain is None or a.chain == chain]
    if not prot or not nuc:
        core = _profile_pairwise([], [], _classify_protein_nucleic)
    else:
        stacking = _ring_stacking(
            prot, AROMATIC_RINGS, nuc, NUCLEIC_BASE_RINGS, STACK_CENTROID_MAX
        )
        core = _profile_pairwise(
            prot, nuc, _classify_protein_nucleic, extra=stacking
        )
    return {
        "mode": "protein-nucleic",
        "protein_chains": sorted({a.chain for a in prot}),
        "nucleic_chains": sorted({a.chain for a in nuc}),
        "interactions": core["interactions"],
        "counts": core["counts"],
        "interaction_total": core["interaction_total"],
        "protein_residues": core["residues_a"],
        "nucleic_residues": core["residues_b"],
        "contact_residue_count": len(core["residues_a"]) + len(core["residues_b"]),
    }
