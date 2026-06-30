"""Antibody CDR loop assignment (Kabat / IMGT numbering schemes).

Annotates the paratope residues from ``interactions.profile_interface`` with
CDR loop labels (H1–H3 for heavy chain, L1–L3 for light chain) based on the
residue-sequence numbers already present in the PDB ATOM records.

Background
----------
Many antibody structures in the PDB are deposited using Kabat positional
numbering (the depositor re-labels residues to match the canonical scheme).
When this convention is followed, standard Kabat windows (31–35, 50–65,
95–102 for VH; 24–34, 50–56, 89–97 for VL) map directly onto ``res_seq``
without re-numbering.  For structures deposited with plain sequential numbering
the CDR windows will not align, so we report a ``numbering_note`` warning in
that case.

Scope: structural interpretation only.  No sequence alignment, no ML
re-numbering (ANARCI), no affinity prediction.  Research-only.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# CDR boundary tables
# ---------------------------------------------------------------------------

# Kabat (KABAT / Wu-Kabat) — the most widely used scheme in the PDB.
_KABAT_H = {
    "H1": (31, 35),
    "H2": (50, 65),
    "H3": (95, 102),
}
_KABAT_L = {
    "L1": (24, 34),
    "L2": (50, 56),
    "L3": (89, 97),
}

# IMGT (the other standard; positions are the same across VH and VL).
_IMGT_H = {
    "H1": (27, 38),
    "H2": (56, 65),
    "H3": (105, 117),
}
_IMGT_L = {
    "L1": (27, 38),
    "L2": (56, 65),
    "L3": (105, 117),
}

# Chothia (structural loop definition — closest to buried cavity geometry).
_CHOTHIA_H = {
    "H1": (26, 32),
    "H2": (52, 56),
    "H3": (95, 102),
}
_CHOTHIA_L = {
    "L1": (24, 34),
    "L2": (50, 56),
    "L3": (89, 97),
}

SCHEMES: dict[str, dict[str, dict]] = {
    "kabat":   {"H": _KABAT_H,   "L": _KABAT_L},
    "imgt":    {"H": _IMGT_H,    "L": _IMGT_L},
    "chothia": {"H": _CHOTHIA_H, "L": _CHOTHIA_L},
}

# Chain-length heuristics (residue counts, not atom counts).
# Fab heavy chain includes VH + CH1 (~215–230 aa).
# Fab light chain includes VL + CL (~210–220 aa, kappa or lambda).
# Isolated variable domain (VH or VL alone): ~110–130 aa.
_FAB_HEAVY = (190, 260)
_FAB_LIGHT = (180, 240)
_VAR_DOMAIN = (100, 140)


def _cdr_for_pos(res_seq: int, table: dict) -> str | None:
    for name, (lo, hi) in table.items():
        if lo <= res_seq <= hi:
            return name
    return None


def _framework_coverage(res_seqs: set, table: dict) -> int:
    """Count how many CDR windows have at least one residue present."""
    return sum(
        1 for lo, hi in table.values()
        if any(lo <= r <= hi for r in res_seqs)
    )


# ---------------------------------------------------------------------------
# Chain-length-based auto-detection
# ---------------------------------------------------------------------------

def detect_vhvl_chains(structure) -> dict:
    """Heuristically identify the heavy (VH) and light (VL) chains.

    Uses chain residue-count ranges typical of Fab fragments or isolated
    variable domains, then scores how many CDR windows are covered by the
    actual residue numbers to prefer chains with Kabat-style numbering.

    Returns a dict with keys: heavy, light, confidence, notes.
    """
    by_chain: dict[str, set] = {}
    for a in structure.protein_atoms:
        by_chain.setdefault(a.chain, set()).add(a.res_seq)
    lengths = {c: len(s) for c, s in by_chain.items()}

    fab_h_cands = [(c, n) for c, n in lengths.items()
                   if _FAB_HEAVY[0] <= n <= _FAB_HEAVY[1]]
    fab_l_cands = [(c, n) for c, n in lengths.items()
                   if _FAB_LIGHT[0] <= n <= _FAB_LIGHT[1]]
    vd_cands    = [(c, n) for c, n in lengths.items()
                   if _VAR_DOMAIN[0] <= n <= _VAR_DOMAIN[1]]

    notes: list[str] = []

    def pick(cands, kabat_table):
        def score(c):
            return _framework_coverage(by_chain[c], kabat_table)
        return max(cands, key=lambda t: (score(t[0]), t[1]))[0]

    heavy = light = None

    if fab_h_cands and fab_l_cands:
        heavy = pick(fab_h_cands, _KABAT_H)
        l_cands = [c for c in fab_l_cands if c[0] != heavy]
        light = pick(l_cands, _KABAT_L) if l_cands else None
        conf = "high"
        notes.append(
            f"Fab-length chains detected (H~{lengths[heavy]}, L~{lengths[light]} residues)"
            if light else f"Fab heavy chain detected (H~{lengths[heavy]} residues)"
        )
    elif fab_h_cands:
        heavy = pick(fab_h_cands, _KABAT_H)
        conf = "medium"
        notes.append(f"Heavy-chain-length chain found (~{lengths[heavy]} residues); no matching light chain")
    elif fab_l_cands:
        light = pick(fab_l_cands, _KABAT_L)
        conf = "medium"
        notes.append(f"Light-chain-length chain found (~{lengths[light]} residues); no matching heavy chain")
    elif len(vd_cands) >= 2:
        heavy, light = vd_cands[0][0], vd_cands[1][0]
        conf = "low"
        notes.append(
            "Two variable-domain-length chains (~110–130 aa) — may be an scFv; "
            "H/L assignment is tentative"
        )
    else:
        conf = "low"
        notes.append("No chains match typical antibody Fab or variable-domain dimensions")

    return {"heavy": heavy, "light": light, "confidence": conf, "notes": notes}


# ---------------------------------------------------------------------------
# Paratope CDR annotation
# ---------------------------------------------------------------------------

def _numbering_looks_kabat(res_seqs: set) -> bool:
    """Check whether residue numbers plausibly come from Kabat renumbering.

    Kabat VH numbers rarely exceed 113; VL numbers rarely exceed 107 (kappa)
    or 110 (lambda).  Structures with sequential numbering from 1 to the total
    chain length will far exceed these bounds for typical Fab chains.
    """
    return bool(res_seqs) and max(res_seqs) <= 130


def annotate_paratope(
    interface_profile: dict,
    heavy: str | None,
    light: str | None,
    scheme: str = "kabat",
) -> dict:
    """Annotate paratope (side-A) residues with CDR loop labels.

    Adds a ``cdr`` key to each residue in ``interface_residues_a`` that falls
    within a CDR window; ``None`` for framework residues.  Also adds:

    * ``cdr_summary``       — {CDR_name: contact_count} for the paratope.
    * ``antibody_chains``   — {heavy, light, scheme}.
    * ``numbering_note``    — warning if the residue numbers look non-Kabat.
    * ``has_cdr_annotation`` — True.

    Returns a shallow copy of the profile dict (original not mutated).
    """
    if scheme not in SCHEMES:
        scheme = "kabat"
    h_table = SCHEMES[scheme]["H"]
    l_table = SCHEMES[scheme]["L"]

    residues_a = [dict(r) for r in (interface_profile.get("interface_residues_a") or [])]
    cdr_summary: dict[str, int] = {}
    all_seqs: set[int] = set()

    for r in residues_a:
        chain = r.get("chain")
        seq = r.get("res_seq", -1)
        all_seqs.add(seq)
        if chain == heavy:
            cdr = _cdr_for_pos(seq, h_table)
        elif chain == light:
            cdr = _cdr_for_pos(seq, l_table)
        else:
            cdr = None
        r["cdr"] = cdr
        if cdr:
            cdr_summary[cdr] = cdr_summary.get(cdr, 0) + r.get("total", 1)

    note = None
    if (heavy or light) and not _numbering_looks_kabat(all_seqs):
        note = (
            "Residue numbers in the paratope chains exceed typical Kabat range "
            f"(max={max(all_seqs) if all_seqs else 'n/a'}).  This structure may "
            "use sequential (non-Kabat) numbering — CDR windows may not align "
            "correctly.  Deposit or renumber with ANARCI for reliable CDR calls."
        )

    profile = dict(interface_profile)
    profile["interface_residues_a"] = residues_a
    profile["cdr_summary"] = cdr_summary
    profile["antibody_chains"] = {"heavy": heavy, "light": light, "scheme": scheme}
    profile["numbering_note"] = note
    profile["has_cdr_annotation"] = True
    return profile
