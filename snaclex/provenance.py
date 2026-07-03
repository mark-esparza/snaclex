"""Method-transparency / provenance blocks for each analysis.

The audit's single "highest-leverage product idea" is a benchmark-first mode:
every algorithm surface should expose its method family, version, key
parameters, a plain-language "what this score means / does not mean", and its
limitations. Docking and screening already emit a methods block (see
``server._methods_block``); this module adds matching blocks for pocket
detection and conservation so *every* analytical result carries provenance.

These blocks are static for a given build (they describe the algorithm, not a
particular structure), which keeps them trivially unit-testable.
"""

from __future__ import annotations

import json
import os

from . import __version__
from . import evolution, pockets

_BENCH_PATH = os.path.join(os.path.dirname(__file__), os.pardir, "benchmark_results.json")

_RESEARCH_DISCLAIMER = (
    "Research-only. A heuristic computed from a single static structure (and, "
    "for conservation, a family alignment). Not an affinity, not a validated "
    "binding-site or functional-site call, not clinical guidance. Validate with "
    "orthogonal evidence."
)


def tool() -> str:
    return f"SnaCleX v{__version__}"


def docking_benchmark(path=None):
    """Return the latest committed docking-benchmark summary, or None.

    Populated by committing the output of ``python -m snaclex.benchmark`` to
    ``benchmark_results.json`` at the repo root, so every docking result can
    link to the method's last measured pose-recovery numbers.
    """
    try:
        with open(path or _BENCH_PATH, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    s = data.get("summary") or {}
    return {
        "last_benchmark_utc": data.get("run_utc"),
        "datasets": data.get("datasets"),
        "n_cases": s.get("n_cases"),
        "top1_success_rate": s.get("top1_success_rate"),
        "median_rmsd_A": s.get("median_rmsd_A"),
        "note": (
            "Self-docking RMSD recovery (nearest-atom proxy); sub-2 Å counts as "
            "the crystallographic pose recovered."
        ),
    }


def pocket_methods() -> dict:
    """Provenance for the LIGSITE pocket-detection + druggability scoring."""
    return {
        "tool": tool(),
        "method": "LIGSITE geometric cavity detection (Hendlich et al., 1997)",
        "method_family": "geometric pocket detection",
        "parameters": {
            "psp_threshold": pockets.PSP_THRESHOLD,
            "scan_range_steps": pockets.SCAN_RANGE,
            "min_pocket_points": pockets.MIN_POCKET_POINTS,
            "fill_radius_A": pockets.FILL,
            "box_margin_A": pockets.MARGIN,
            "max_pockets_reported": 6,
        },
        "scoring": (
            "Druggability index (0-100) = 0.30·volume + 0.30·enclosure + "
            "0.30·hydrophobicity − 0.15·polarity (+0.15 offset), clamped to "
            "[0,100]. A transparent rule-of-thumb proxy for SiteMap/fpocket-style "
            "druggability — NOT a trained model."
        ),
        "interpretation": (
            "Tiers: 'pocket' = geometric cavity only; 'ligandable' = large and "
            "enclosed enough to bind a small molecule; 'druggable' = also has "
            "favourable (hydrophobic, enclosed) chemistry. A higher score means a "
            "more promising cavity in relative terms; it does NOT confirm a real "
            "binding site or predict affinity."
        ),
        "limitations": [
            "Pure geometry on a single conformation — no induced fit or "
            "cryptic-pocket dynamics.",
            "Grid-discretized: very shallow or very small cavities may be missed "
            "or merged.",
            "Druggability is a heuristic combination of volume, enclosure and "
            "residue chemistry, not a benchmarked classifier.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }


def antibody_methods() -> dict:
    """Provenance for antibody chain typing, CDR delimitation + liability scan."""
    from . import antibody
    return {
        "tool": tool(),
        "method": "Germline-framework sequence alignment (chain typing + CDR "
                  "delimitation) with regex developability-liability scan",
        "method_family": "antibody sequence analysis",
        "parameters": {
            "numbering_scheme": "IMGT (CDR delimitation)",
            "typing": "Needleman-Wunsch alignment of each chain's N-terminal "
                      "window to VH / V-kappa / V-lambda / CH1 / CL germline "
                      "framework consensus references",
            "variable_framework_identity_min": antibody._VAR_FR_MIN,
            "anchor_check": "two invariant intradomain cysteines must align",
            "liability_motifs": {
                "N-glycosylation": "N-X-S/T (X≠Pro)",
                "deamidation": "N-G / N-S",
                "isomerization": "D-G",
                "unpaired_cysteine": "Cys beyond the two canonical anchors / odd count",
            },
        },
        "interpretation": (
            "Chain typing and CDR ranges are a fast germline-alignment heuristic "
            "(not ANARCI/IMGT-HMM numbering). Liabilities are sequence motifs; a "
            "motif in a CDR is higher-concern than one in framework. These are "
            "developability flags, NOT a stability or immunogenicity prediction."
        ),
        "limitations": [
            "Framework-alignment typing can mis-call unusual germlines, "
            "single-domain (VHH) or engineered scaffolds; CDR boundaries are "
            "approximate, especially CDR-H3.",
            "Liability motifs are presence/absence only — no structural exposure, "
            "pH, or formulation context is considered.",
            "Interface buried-surface-area is a coarse Shrake-Rupley estimate "
            "(heavy atoms, no hydrogens), not a rigorous SASA.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }


def evolution_methods() -> dict:
    """Provenance for the Pfam-alignment conservation analysis."""
    return {
        "tool": tool(),
        "method": "Pfam/InterPro family-alignment conservation",
        "method_family": "evolutionary conservation",
        "parameters": {
            "alignment_source": "Pfam family alignment via InterPro / EMBL-EBI",
            "max_alignment_sequences": evolution.MAX_ALIGN_SEQS,
            "conservation_metric": "per-column Shannon entropy, 1 − H/log2(20)",
            "mapping": "Needleman-Wunsch global alignment of the structure "
                       "sequence onto alignment columns",
        },
        "interpretation": (
            "Per-residue conservation runs 0 (variable) to 1 (invariant) from "
            "the family alignment column entropy; a pocket's conservation is the "
            "mean over its lining residues. High conservation suggests functional "
            "or structural importance — it does NOT by itself mean a residue is "
            "part of a ligand-binding site."
        ),
        "limitations": [
            "Depends on the protein mapping to a Pfam family with a usable "
            "alignment; unavailable otherwise.",
            "Shallow or skewed alignments give noisy per-column entropy.",
            "Conservation reflects evolutionary pressure broadly (fold, catalysis, "
            "interfaces), not binding specifically.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }
