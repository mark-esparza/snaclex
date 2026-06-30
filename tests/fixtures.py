"""Builders for small, deterministic, offline test structures.

These construct ``snaclex.pdbparse`` dataclasses directly so the compute modules
(interactions, docking, pockets) can be exercised without any network access or
real PDB files.
"""

from __future__ import annotations

from snaclex.pdbparse import Atom, Component, Structure


def atom(
    element,
    x,
    y,
    z,
    *,
    name=None,
    res_name="ALA",
    chain="A",
    res_seq=1,
    hetero=False,
    serial=1,
    icode="",
):
    """Build a single Atom with sensible defaults."""
    return Atom(
        serial=serial,
        name=name or element,
        res_name=res_name,
        chain=chain,
        res_seq=res_seq,
        icode=icode,
        x=x,
        y=y,
        z=z,
        element=element,
        is_hetero=hetero,
    )


def component(res_name, atoms, *, chain="X", res_seq=900, icode=""):
    """Build a hetero Component from a list of atoms."""
    return Component(res_name, chain, res_seq, icode, list(atoms))


_MMCIF_TAGS = [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code",
    "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
    "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id",
    "pdbx_PDB_model_num",
]


def mmcif_text(rows):
    """Build a minimal mmCIF `_atom_site` loop from row dicts.

    Each row may set: group, element, name, comp, chain, seq, x, y, z, alt,
    icode, model. Sensible defaults fill the rest.
    """
    out = ["data_TEST", "#", "loop_"]
    out += [f"_atom_site.{t}" for t in _MMCIF_TAGS]
    for i, r in enumerate(rows, start=1):
        name = r.get("name", "CA")
        comp = r.get("comp", "LEU")
        chain = r.get("chain", "A")
        seq = r.get("seq", i)
        vals = [
            r.get("group", "ATOM"), i, r.get("element", "C"), name,
            r.get("alt", "."), comp, chain, seq, r.get("icode", "?"),
            r.get("x", 0.0), r.get("y", 0.0), r.get("z", 0.0),
            r.get("occ", 1.0), r.get("b", 20.0),
            seq, comp, chain, name, r.get("model", 1),
        ]
        out.append(" ".join(str(v) for v in vals))
    out.append("#")
    return "\n".join(out) + "\n"


def structure(protein=None, components=None, nucleic=None):
    """Assemble a Structure from protein atoms, hetero components, nucleic atoms."""
    protein = list(protein or [])
    components = list(components or [])
    nucleic = list(nucleic or [])
    all_atoms = protein + nucleic + [a for c in components for a in c.atoms]
    chains = sorted({a.chain for a in protein})
    nucleic_chains = sorted({a.chain for a in nucleic})
    return Structure(
        atoms=all_atoms,
        protein_atoms=protein,
        components=components,
        chains=chains,
        nucleic_atoms=nucleic,
        nucleic_chains=nucleic_chains,
    )
