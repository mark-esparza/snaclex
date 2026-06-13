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


def structure(protein=None, components=None):
    """Assemble a Structure from protein atoms and hetero components."""
    protein = list(protein or [])
    components = list(components or [])
    all_atoms = protein + [a for c in components for a in c.atoms]
    chains = sorted({a.chain for a in protein})
    return Structure(
        atoms=all_atoms,
        protein_atoms=protein,
        components=components,
        chains=chains,
    )
