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


def interface_methods() -> dict:
    """Provenance for protein-protein / protein-nucleic interface profiling."""
    return {
        "tool": tool(),
        "method": "Geometric heavy-atom interface contact classification",
        "method_family": "interaction profiling",
        "parameters": {
            "hydrogen_bond_max_A": 3.6,
            "salt_bridge_max_A": 4.0,
            "hydrophobic_range_A": "2.8-4.0",
            "ring_stacking_centroid_max_A": 6.0,
            "engine": "shared spatial grid (same as ligand interaction profiler)",
        },
        "interpretation": (
            "Contacts between two chain groups (antibody-antigen, HLA-peptide) or "
            "between protein and nucleic-acid chains are classified into H-bonds, "
            "salt bridges (incl. basic side chain .. DNA/RNA phosphate), hydrophobic "
            "contacts and ring stacking, then summarized per residue to give an "
            "epitope/paratope (or protein/nucleic) contact map. It locates the "
            "interface; it does NOT compute binding energy or affinity."
        ),
        "limitations": [
            "Heavy-atom geometry on one static conformation; no explicit hydrogens, "
            "no induced fit, no energetics.",
            "Antibody CDR numbering is not applied (interface residues only).",
            "Chain grouping is user/heuristic-driven; mis-grouped chains mislead.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }


def variant_methods() -> dict:
    """Provenance for variant -> structure-residue mapping."""
    return {
        "tool": tool(),
        "method": "UniProt-to-structure variant localization",
        "method_family": "variant mapping",
        "parameters": {
            "alignment": "Needleman-Wunsch global alignment (UniProt canonical "
                         "sequence vs structure sequence)",
            "input": "protein position (R273H / p.Arg273His / bare position); "
                     "genomic chr:pos ref>alt via Ensembl VEP when enabled",
            "population_frequency": "TOPMed/BRAVO aggregate allele frequency "
                                    "(public; genomic path only)",
        },
        "interpretation": (
            "Each variant is placed on a concrete (chain, residue) in the 3D model "
            "and cross-referenced with pocket / interface / conservation context. "
            "This is a structural localization, NOT a pathogenicity or clinical call."
        ),
        "limitations": [
            "Only positions modeled in the structure can be mapped; unmodeled "
            "regions and isoform differences are flagged, not mapped.",
            "BRAVO frequencies are public aggregates only — never individual "
            "genotypes or controlled-access data.",
            "WT/structure mismatches (engineered or alternate isoforms) are flagged "
            "but not corrected.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }


def hla_methods() -> dict:
    """Provenance for HLA peptide-binding-groove analysis."""
    return {
        "tool": tool(),
        "method": "HLA/MHC groove + anchor-pocket structural analysis",
        "method_family": "immunogenomics (structural)",
        "parameters": {
            "detection": "chain composition (heavy + beta-2-microglobulin + short "
                         "peptide for class I) plus title/UniProt signals",
            "groove": "HLA residues contacting the bound peptide (interface profiler)",
            "anchors": "class-I P2 (B-pocket) and C-terminal PΩ (F-pocket)",
            "variant_classes": ["anchor-pocket", "groove-lining", "peripheral"],
        },
        "interpretation": (
            "Identifies the peptide-binding groove and anchor pockets and classifies "
            "variants by where they fall, giving a transparent 'likely alters the "
            "peptide repertoire vs peripheral' read-out. This is STRUCTURAL "
            "INTERPRETATION ONLY — not peptide-binding-affinity or immunogenicity "
            "prediction (use NetMHCpan-style tools for that) and not HLA typing."
        ),
        "limitations": [
            "Requires a modeled peptide; apo structures give detection only.",
            "Anchor assignment is class-I-centric (P2 / C-terminus); class-II "
            "anchors are not sub-classified.",
            "Allele-specific repertoire effects are not quantified.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }


def esm_methods() -> dict:
    """Provenance for the optional EvolutionaryScale ESM model layer."""
    return {
        "tool": tool(),
        "method": "EvolutionaryScale ESM (Forge): ESM3 fold + ESMC variant scoring",
        "method_family": "protein language model (external, optional)",
        "parameters": {
            "fold_model": "esm3-open (structure from sequence)",
            "score_model": "esmc-600m (masked-marginal log-likelihood ratio)",
            "provider": "EvolutionaryScale Forge (requires ESM_API_KEY)",
        },
        "interpretation": (
            "Predicted structures let the structural pipeline run when no "
            "experimental model exists; variant log-likelihood ratios provide a "
            "sequence-model 'predictability' signal shown ALONGSIDE the structural "
            "overlay. Both are model PREDICTIONS, labelled as such."
        ),
        "limitations": [
            "Off unless a Forge token is configured; the default build is "
            "structural-only and unaffected.",
            "Predicted structures carry per-region uncertainty and are not "
            "experimental coordinates.",
            "Likelihood ratios are not binding-affinity or pathogenicity values.",
        ],
        "disclaimer": _RESEARCH_DISCLAIMER,
    }
