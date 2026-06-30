"""Heuristic antibody CDR (complementarity-determining region) detection.

Identifies candidate immunoglobulin variable-domain (Fv) chains by sequence
length, then locates the three CDR loops relative to the framework anchors
that flank them: the Cys ~22/23 (Kabat) ending FR1, the conserved Trp opening
FR2, the Cys ~92 (Kabat) ending FR3, and the conserved "W-G-x-G" (heavy) /
"F-G-x-G" (light) motif opening FR4.

This is a lightweight, alignment-free heuristic -- not a substitute for a
full germline-aligned numbering scheme (ANARCI/IMGT/Chothia). CDR3, anchored
on the well-conserved Cys...[W/F]Gx G framework motif, is the most reliable;
CDR1/CDR2 boundaries are approximate fixed offsets from the Trp anchor.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .pdbparse import Structure

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}

_MIN_V_DOMAIN = 95
_MAX_V_DOMAIN = 140

# Conserved FR3 -> FR4 framework motif: heavy "WGQG"-like, light "FGxG"-like.
_J_MOTIF = re.compile(r"WG.G|FG.G")


@dataclass
class CDRLoop:
    name: str
    chain: str
    start_res_seq: int
    end_res_seq: int
    sequence: str
    confidence: str  # "high" | "approximate"


def _chain_sequence(structure: Structure, chain: str) -> tuple[list[int], str]:
    """CA-ordered (residue numbers, one-letter sequence) for one protein chain."""
    residues: dict[int, str] = {}
    for a in structure.protein_atoms:
        if a.chain != chain or a.name != "CA":
            continue
        residues[a.res_seq] = AA3TO1.get(a.res_name, "X")
    res_seqs = sorted(residues)
    seq = "".join(residues[r] for r in res_seqs)
    return res_seqs, seq


def _find_cdrs_in_chain(chain: str, res_seqs: list[int], seq: str) -> list[CDRLoop]:
    if not (_MIN_V_DOMAIN <= len(seq) <= _MAX_V_DOMAIN):
        return []

    # FR1 conserved Cys (~Kabat 22/23), expected within the first third.
    c1 = next((i for i in range(15, min(35, len(seq))) if seq[i] == "C"), None)
    if c1 is None:
        return []

    # FR2 conserved Trp, a few residues after c1.
    w1 = next(
        (i for i in range(c1 + 5, min(c1 + 25, len(seq))) if seq[i] == "W"), None
    )
    if w1 is None:
        return []

    # FR3 conserved Cys (~Kabat 92), 55-85 residues after c1.
    c2 = next(
        (i for i in range(c1 + 55, min(c1 + 85, len(seq))) if seq[i] == "C"), None
    )
    if c2 is None:
        return []

    # FR4 anchor motif, after c2. CDR3 itself can contain a spurious match of
    # this pattern, so take the last match in the chain (FR4 sits near the
    # domain's C-terminus) rather than the first.
    matches = list(_J_MOTIF.finditer(seq, c2 + 3))
    if not matches:
        return []
    j_start = matches[-1].start()
    if j_start <= c2 + 2:
        return []

    def loop(name: str, start: int, end: int, confidence: str) -> CDRLoop | None:
        start = max(0, start)
        end = min(len(seq) - 1, end)
        if end < start:
            return None
        return CDRLoop(
            name=name,
            chain=chain,
            start_res_seq=res_seqs[start],
            end_res_seq=res_seqs[end],
            sequence=seq[start : end + 1],
            confidence=confidence,
        )

    loops = [
        loop("CDR1", c1 + 4, w1 - 3, "approximate"),
        loop("CDR2", w1 + 14, w1 + 24, "approximate"),
        loop("CDR3", c2 + 3, j_start - 1, "high"),
    ]
    return [item for item in loops if item is not None]


def detect_cdrs_json(structure: Structure) -> list[dict]:
    """``detect_cdrs`` results as JSON-serializable dicts."""
    return [
        {
            "name": loop.name,
            "chain": loop.chain,
            "start_res_seq": loop.start_res_seq,
            "end_res_seq": loop.end_res_seq,
            "sequence": loop.sequence,
            "confidence": loop.confidence,
        }
        for loop in detect_cdrs(structure)
    ]


def detect_cdrs(structure: Structure) -> list[CDRLoop]:
    """Detect candidate antibody CDR loops across all protein chains.

    Scans each chain for an immunoglobulin-variable-domain-sized sequence
    bracketed by the conserved Cys...Cys...[W/F]Gx G framework anchors and
    reports the CDR1/2/3 loops relative to those anchors. Chains that don't
    match the expected size/anchor pattern (i.e. aren't antibody Fv domains)
    are skipped, so this is safe to call on any structure.
    """
    out: list[CDRLoop] = []
    for chain in structure.chains:
        res_seqs, seq = _chain_sequence(structure, chain)
        out.extend(_find_cdrs_in_chain(chain, res_seqs, seq))
    return out
