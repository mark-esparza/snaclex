"""Sequence-only protein analysis (Stage 3) — pure Python, no network.

Every metric here is *calculated* from the amino-acid sequence alone, so it is
available for any protein regardless of whether a structure exists. Results are
labelled ``evidence_category = "calculated"`` by the record builder.

Nothing in this module fetches anything; it is deterministic and fully unit
testable offline, matching the existing pure-compute modules (pockets, docking,
interactions). The methods are transparent textbook calculations, not trained
models:

* Molecular weight  — sum of average residue masses + one water.
* Theoretical pI    — EMBOSS pKa set, net-charge root found by bisection.
* Net charge at pH  — Henderson-Hasselbalch over ionizable groups.
* GRAVY / hydropathy — Kyte & Doolittle (1982) hydropathy index.
* Low-complexity    — windowed normalized Shannon entropy (a SEG-style proxy).

These are screening heuristics for interpretation, not validated biophysical
measurements.
"""

from __future__ import annotations

import math

# Average isotopic masses of amino-acid *residues* (Da); water added once.
_RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
_WATER = 18.01524

# Kyte & Doolittle (1982) hydropathy index.
_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# EMBOSS pKa set. Documented so the assumption is inspectable.
_PKA_NTERM = 8.6
_PKA_CTERM = 3.6
_PKA_POS = {"K": 10.8, "R": 12.5, "H": 6.5}        # protonated = positive
_PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1}  # deprotonated = negative

_STANDARD = set(_RESIDUE_MASS)
# Non-standard / ambiguity codes we tolerate but exclude from mass/charge.
_AMBIGUOUS = set("XBZJUO*")

PKA_SET_NAME = "EMBOSS"


class SequenceError(ValueError):
    """Raised when the input is not a usable amino-acid sequence."""


def clean_sequence(seq: str) -> str:
    """Uppercase, strip whitespace/digits/gaps; keep only letter codes.

    Raises ``SequenceError`` if nothing usable remains or the residue alphabet
    is mostly invalid (guards against a nucleotide or free-text paste).
    """
    if not seq:
        raise SequenceError("empty sequence")
    cleaned = "".join(ch for ch in seq.upper() if ch.isalpha() or ch == "*")
    if not cleaned:
        raise SequenceError("no amino-acid letters in input")
    valid = sum(1 for ch in cleaned if ch in _STANDARD or ch in _AMBIGUOUS)
    if valid / len(cleaned) < 0.9:
        raise SequenceError("input does not look like an amino-acid sequence")
    return cleaned


def composition(seq: str) -> dict:
    """Return per-residue counts and fractions plus non-standard tallies."""
    counts: dict[str, int] = {}
    for ch in seq:
        counts[ch] = counts.get(ch, 0) + 1
    n = len(seq)
    standard = {a: counts.get(a, 0) for a in sorted(_STANDARD)}
    fractions = {a: (standard[a] / n if n else 0.0) for a in standard}
    non_standard = {a: c for a, c in counts.items() if a not in _STANDARD}
    return {
        "length": n,
        "counts": standard,
        "fractions": fractions,
        "non_standard": non_standard,
    }


def molecular_weight(seq: str) -> float:
    """Average molecular weight in Da (standard residues + one water).

    Non-standard residues are skipped and noted by ``composition``; the value is
    therefore a lower-bound approximation when ambiguity codes are present.
    """
    mass = sum(_RESIDUE_MASS[a] for a in seq if a in _RESIDUE_MASS)
    return round(mass + _WATER, 2) if mass else 0.0


def net_charge(seq: str, ph: float) -> float:
    """Net charge at a given pH (Henderson-Hasselbalch over ionizable groups)."""
    pos = 1.0 / (1.0 + 10 ** (ph - _PKA_NTERM))           # N-terminus
    neg = 1.0 / (1.0 + 10 ** (_PKA_CTERM - ph))           # C-terminus
    for a, pka in _PKA_POS.items():
        pos += seq.count(a) * (1.0 / (1.0 + 10 ** (ph - pka)))
    for a, pka in _PKA_NEG.items():
        neg += seq.count(a) * (1.0 / (1.0 + 10 ** (pka - ph)))
    return pos - neg


def isoelectric_point(seq: str) -> float:
    """Theoretical pI: the pH where net charge crosses zero (bisection)."""
    lo, hi = 0.0, 14.0
    # Net charge is monotonically decreasing in pH, so bisection is exact enough.
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def gravy(seq: str) -> float | None:
    """Grand average of hydropathy (mean KD index over standard residues)."""
    vals = [_KD[a] for a in seq if a in _KD]
    return round(sum(vals) / len(vals), 3) if vals else None


def hydropathy_profile(seq: str, window: int = 9) -> list[float]:
    """Sliding-window mean KD hydropathy, one value per window centre.

    Returns an empty list if the sequence is shorter than the window.
    """
    if window < 1 or len(seq) < window:
        return []
    vals = [_KD.get(a, 0.0) for a in seq]
    out: list[float] = []
    running = sum(vals[:window])
    out.append(round(running / window, 3))
    for i in range(window, len(vals)):
        running += vals[i] - vals[i - window]
        out.append(round(running / window, 3))
    return out


def _window_entropy(counts: dict[str, int], size: int) -> float:
    """Normalized Shannon entropy (0..1) of a window's residue distribution."""
    h = 0.0
    for c in counts.values():
        if c:
            p = c / size
            h -= p * math.log2(p)
    return h / math.log2(20)  # normalize against the 20-letter alphabet


def low_complexity_regions(seq: str, window: int = 20,
                           threshold: float = 0.5) -> list[dict]:
    """Flag windows whose normalized composition entropy falls below threshold.

    A transparent SEG-style proxy: low entropy => repetitive / biased composition
    (e.g. poly-Q, poly-A). 1-based inclusive ``start``/``end``. Heuristic, not a
    substitute for SEG/fLPS.
    """
    n = len(seq)
    if n < window:
        return []
    flagged: list[tuple[int, int]] = []
    for i in range(0, n - window + 1):
        w = seq[i:i + window]
        counts: dict[str, int] = {}
        for ch in w:
            counts[ch] = counts.get(ch, 0) + 1
        if _window_entropy(counts, window) < threshold:
            flagged.append((i + 1, i + window))
    # Merge overlapping windows into contiguous regions.
    merged: list[dict] = []
    for start, end in flagged:
        if merged and start <= merged[-1]["end"] + 1:
            merged[-1]["end"] = max(merged[-1]["end"], end)
        else:
            merged.append({"start": start, "end": end})
    return merged


def analyze(seq: str, ph: float = 7.0, *, profile_points: int = 200) -> dict:
    """Full Stage-3 sequence analysis for an amino-acid sequence.

    ``ph`` sets the pH for the estimated net charge. ``profile_points`` caps the
    returned hydropathy profile (downsampled) so long proteins stay compact.
    All values are *calculated*, labelled as such by the caller.
    """
    seq = clean_sequence(seq)
    comp = composition(seq)
    full_profile = hydropathy_profile(seq)
    if profile_points and len(full_profile) > profile_points:
        step = math.ceil(len(full_profile) / profile_points)
        profile = full_profile[::step]
    else:
        profile = full_profile
    return {
        "length": len(seq),
        "molecular_weight_Da": molecular_weight(seq),
        "isoelectric_point": isoelectric_point(seq),
        "charge_at_ph": {"ph": ph, "net_charge": round(net_charge(seq, ph), 2)},
        "gravy": gravy(seq),
        "aromaticity": round(
            sum(seq.count(a) for a in "FWY") / len(seq), 4) if seq else None,
        "composition": comp,
        "hydropathy_profile": {"window": 9, "downsampled": profile,
                               "n_points_full": len(full_profile)},
        "low_complexity_regions": low_complexity_regions(seq),
        "method": {
            "molecular_weight": "sum of average residue masses + water",
            "isoelectric_point": f"{PKA_SET_NAME} pKa set, bisection on net charge",
            "gravy": "Kyte & Doolittle (1982) hydropathy index",
            "low_complexity": "windowed normalized Shannon entropy (SEG-style proxy)",
            "note": "Calculated screening heuristics — not measured biophysical values.",
        },
    }
