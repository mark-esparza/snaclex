"""Antibody detection, CDR delimitation and developability-liability scan.

Dependency-free. Runs after a structure loads (cheap sequence work) to answer:
is this an antibody, which chains are the variable domains, where are the CDR
loops, and does the sequence carry known developability-liability motifs.

Approach (deliberately lightweight, per the brief):
  * Each protein chain's N-terminal sequence is Needleman-Wunsch aligned (reusing
    ``evolution._nw_align``) to a small set of germline **framework consensus**
    references — VH, V-kappa, V-lambda and the constant CH1 / CL domains.
  * A chain is typed VH / VL only when its framework identity is high enough AND
    the two invariant immunoglobulin-domain cysteines line up. Requiring the
    conserved disulfide cysteines (not just fold similarity) keeps this from
    firing on every Ig-superfamily protein (e.g. beta-2-microglobulin, CD-fold
    antigens).
  * CDR spans are read off the reference's annotated CDR columns; CDR3's variable
    length is captured by walking the alignment and attaching every inserted
    residue to the last framework segment seen, then snapped to the J-region
    W/F-G-x-G motif.

Numbering scheme reported is IMGT (default). Research-only; this is a fast
heuristic, not ANARCI/IMGT-HMM numbering.
"""

from __future__ import annotations

import re

from .evolution import AA3TO1, _nw_align

# Only the N-terminal window of each chain is used for variable-domain typing
# (the V domain is N-terminal, ~110-130 aa); this bounds the NW cost per chain.
_VAR_WINDOW = 160

# Framework-identity thresholds. Framework identity is matches / (framework
# columns of the reference), so unaligned framework columns count against a
# candidate — non-antibodies score low.
# Real V-domains score ~0.6-0.97 framework identity to their own germline
# reference; unrelated folds (even other Ig-superfamily / cysteine-rich proteins)
# stay <=~0.4. 0.45 separates them with margin.
_VAR_FR_MIN = 0.45
_CONST_FR_MIN = 0.45

# ---------------------------------------------------------------------------
# Germline framework consensus references.
#
# Each reference is a list of (segment, subsequence). CDR subsequences are the
# germline placeholders used only for alignment anchoring; their columns are NOT
# scored for identity (they are hypervariable). The two invariant cysteines live
# in FR1 (…S-C) and FR3 (…YY-C) and are the only cysteines in the framework
# strings, so we locate the anchor columns simply as the framework 'C' columns.
# ---------------------------------------------------------------------------

_VH = [
    ("FR1", "EVQLVESGGGLVQPGGSLRLSCAAS"),
    ("CDR1", "GFTFSSYAMS"),
    ("FR2", "WVRQAPGKGLEWVS"),
    ("CDR2", "AISGSGGSTY"),
    # FR3 ends at the second conserved cysteine (IMGT 104); the germline "AK"
    # that follows is IMGT 105-106 and belongs to CDR3, so keep it there.
    ("FR3", "YADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYC"),
    ("CDR3", "AKDRGGYYFDY"),
    ("FR4", "WGQGTLVTVSS"),
]

_VK = [
    ("FR1", "DIQMTQSPSSLSASVGDRVTITC"),
    ("CDR1", "RASQSISSYLN"),
    ("FR2", "WYQQKPGKAPKLLIY"),
    ("CDR2", "AASSLQS"),
    ("FR3", "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYC"),
    ("CDR3", "QQSYSTPPT"),
    ("FR4", "FGQGTKVEIK"),
]

_VL = [
    ("FR1", "QSVLTQPPSVSGAPGQRVTISC"),
    ("CDR1", "TGSSSNIGAGYDVH"),
    ("FR2", "WYQQLPGTAPKLLIY"),
    ("CDR2", "GNSNRPS"),
    ("FR3", "GVPDRFSGSKSGTSASLAITGLQAEDEADYYC"),
    ("CDR3", "QSYDSSLSGSV"),
    ("FR4", "FGGGTKLTVL"),
]

# Constant-domain consensus (single Ig domain each) for typing chains that carry
# no variable domain (pure Fc / CH / CL fragments). No CDRs — all framework.
_CH1 = [("FR", "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV")]
_CL = [("FR", "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGE")]

_REFS_VAR = {"VH": _VH, "VK": _VK, "VL": _VL}
_REFS_CONST = {"heavy_constant": _CH1, "light_constant": _CL}

# Numbering-scheme-agnostic label for what the "VL" bucket really means.
_TYPE_LABEL = {
    "VH": "heavy variable (VH)",
    "VL": "light variable (VL)",
    "heavy_constant": "heavy constant",
    "light_constant": "light constant",
    "other": "other (candidate antigen)",
}

# Developability-liability motifs (sequence-only). Patterns match on the
# one-letter sequence; the reported position is the first residue of the motif.
_LIABILITIES = [
    ("N-glycosylation", re.compile(r"N[^P][ST]"), "N-X-S/T sequon (X≠Pro); potential N-linked glycosylation site."),
    ("Deamidation", re.compile(r"N[GS]"), "Asn-Gly/Ser; hotspot for asparagine deamidation."),
    ("Isomerization", re.compile(r"DG"), "Asp-Gly; hotspot for aspartate isomerization."),
]


def _flatten(ref):
    """Return (ref_seq, seg_per_column) for a segmented reference."""
    seq = []
    segs = []
    for seg, sub in ref:
        for ch in sub:
            seq.append(ch)
            segs.append(seg)
    return "".join(seq), segs


def _chain_sequences(structure):
    """Per chain: (one-letter seq, [res_seq]) in ascending residue order."""
    by_chain: dict[str, dict] = {}
    for a in structure.protein_atoms:
        d = by_chain.setdefault(a.chain, {})
        if a.res_seq not in d:
            d[a.res_seq] = a.res_name
    out = {}
    for chain, residues in by_chain.items():
        items = sorted(residues.items())
        seq = "".join(AA3TO1.get(rn, "X") for _, rn in items)
        out[chain] = (seq, [rs for rs, _ in items])
    return out


def _align_to_ref(seq, ref):
    """Align ``seq`` to a segmented reference.

    Returns a dict with framework identity, whether the anchor cysteines matched,
    and — for variable references — the residue-index spans of each CDR segment.
    Indices are positions into ``seq`` (0-based).
    """
    ref_seq, segs = _flatten(ref)
    fr_cols = sum(1 for s in segs if s.startswith("FR"))
    anchor_cols = [j for j, (c, s) in enumerate(zip(ref_seq, segs))
                   if c == "C" and s.startswith("FR")]

    pairs = _nw_align(seq, ref_seq)

    fr_match = 0
    anchor_hit = 0
    seg_positions: dict[str, list] = {}
    current_seg = None
    for ti, rj in pairs:
        if rj is not None:
            current_seg = segs[rj]
            if segs[rj].startswith("FR") and ti is not None and seq[ti] == ref_seq[rj]:
                fr_match += 1
            if rj in anchor_cols and ti is not None and seq[ti] == "C":
                anchor_hit += 1
        if ti is not None and current_seg is not None:
            seg_positions.setdefault(current_seg, []).append(ti)

    identity = fr_match / fr_cols if fr_cols else 0.0
    return {
        "identity": identity,
        "anchor_hits": anchor_hit,
        "anchor_total": len(anchor_cols),
        "seg_positions": seg_positions,
    }


def _cdr_spans(aln, seq, res_seqs, scheme):
    """Turn a variable-domain alignment into CDR residue spans for one chain."""
    cdrs = {}
    for name in ("CDR1", "CDR2", "CDR3"):
        positions = sorted(aln["seg_positions"].get(name, []))
        if not positions:
            continue
        # Contiguous residue indices only (guard against stray insertions).
        lo, hi = positions[0], positions[-1]
        idxs = list(range(lo, hi + 1))
        cdrs[name] = {
            "res_seqs": [res_seqs[i] for i in idxs],
            "seq": "".join(seq[i] for i in idxs),
            "range": [res_seqs[lo], res_seqs[hi]],
        }
    return cdrs


def _type_chain(seq, res_seqs, scheme):
    """Classify one chain and, if variable, return its CDR spans."""
    window = seq[:_VAR_WINDOW]

    best = None
    for kind, ref in _REFS_VAR.items():
        aln = _align_to_ref(window, ref)
        if best is None or aln["identity"] > best[1]["identity"]:
            best = (kind, aln)

    kind, aln = best
    if aln["identity"] >= _VAR_FR_MIN and aln["anchor_hits"] >= 1:
        vtype = "VH" if kind == "VH" else "VL"
        cdrs = _cdr_spans(aln, window, res_seqs, scheme)
        return {
            "type": vtype,
            "identity": round(aln["identity"], 3),
            "anchor_cys": f"{aln['anchor_hits']}/{aln['anchor_total']}",
            "cdrs": cdrs,
            "germline_ref": kind,
        }

    # No variable domain — try the constant references.
    cbest = None
    for kind, ref in _REFS_CONST.items():
        aln_c = _align_to_ref(seq, ref)
        if cbest is None or aln_c["identity"] > cbest[1]["identity"]:
            cbest = (kind, aln_c)
    if cbest and cbest[1]["identity"] >= _CONST_FR_MIN:
        return {"type": cbest[0], "identity": round(cbest[1]["identity"], 3), "cdrs": {}}

    return {"type": "other", "identity": round(aln["identity"], 3), "cdrs": {}}


def _scan_liabilities(seq, res_seqs, cdr_res_set):
    """Regex liability scan over one chain's sequence."""
    hits = []
    for label, pattern, note in _LIABILITIES:
        for m in pattern.finditer(seq):
            i = m.start()
            pos = res_seqs[i]
            in_cdr = pos in cdr_res_set or (
                i + 1 < len(res_seqs) and res_seqs[i + 1] in cdr_res_set
            )
            hits.append({
                "type": label,
                "motif": m.group(0),
                "res_seq": pos,
                "in_cdr": in_cdr,
                "note": note,
            })

    # Unpaired-cysteine check: variable domains carry the two canonical anchor
    # cysteines; extra ones (or an odd total) suggest a free thiol.
    cys_positions = [res_seqs[i] for i, c in enumerate(seq) if c == "C"]
    if len(cys_positions) > 2:
        for pos in cys_positions[2:]:
            hits.append({
                "type": "Unpaired cysteine",
                "motif": "C",
                "res_seq": pos,
                "in_cdr": pos in cdr_res_set,
                "note": (
                    "Cysteine beyond the two canonical intradomain-disulfide "
                    "cysteines; candidate free thiol (aggregation / heterogeneity)."
                ),
            })
    elif len(cys_positions) % 2 == 1:
        pos = cys_positions[-1]
        hits.append({
            "type": "Unpaired cysteine",
            "motif": "C",
            "res_seq": pos,
            "in_cdr": pos in cdr_res_set,
            "note": "Odd number of cysteines in the variable domain; one is unpaired.",
        })
    return hits


def analyze(structure, scheme: str = "imgt") -> dict:
    """Detect antibody chains, CDRs and liabilities for a loaded structure.

    Always returns a dict; ``is_antibody`` is False when no variable domain is
    found. Cheap enough to run inline on structure load.
    """
    seqs = _chain_sequences(structure)
    chains = []
    liabilities = []
    cdr_residues = []  # flat, for 3D coloring
    any_variable = False
    any_other = False

    for chain in sorted(seqs):
        seq, res_seqs = seqs[chain]
        if len(seq) < 20:
            chains.append({"chain": chain, "type": "other",
                           "label": _TYPE_LABEL["other"], "length": len(seq)})
            any_other = True
            continue
        typed = _type_chain(seq, res_seqs, scheme)
        vtype = typed["type"]
        entry = {
            "chain": chain,
            "type": vtype,
            "label": _TYPE_LABEL[vtype],
            "length": len(seq),
            "framework_identity": typed.get("identity"),
        }
        if vtype in ("VH", "VL"):
            any_variable = True
            entry["anchor_cys"] = typed.get("anchor_cys")
            cdrs = typed.get("cdrs", {})
            entry["cdrs"] = {
                name: {"range": d["range"], "seq": d["seq"], "length": len(d["res_seqs"])}
                for name, d in cdrs.items()
            }
            cdr_res_set = set()
            for name, d in cdrs.items():
                for rs in d["res_seqs"]:
                    cdr_res_set.add(rs)
                    cdr_residues.append({"chain": chain, "res_seq": rs, "cdr": name})
            liabilities.extend(
                {**h, "chain": chain} for h in _scan_liabilities(seq, res_seqs, cdr_res_set)
            )
        elif vtype == "other":
            any_other = True
        chains.append(entry)

    liabilities.sort(key=lambda h: (not h["in_cdr"], h["chain"], h["res_seq"]))

    return {
        "is_antibody": any_variable,
        "scheme": scheme.upper(),
        "chains": chains,
        "cdr_residues": cdr_residues,
        "liabilities": liabilities,
        "liability_count": len(liabilities),
        # A paratope/epitope interface is offerable when an antibody variable
        # chain and a non-antibody ("other") chain coexist in the same file.
        "interface_available": any_variable and any_other,
    }
