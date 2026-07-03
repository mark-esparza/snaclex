"""HLA / MHC module: analyze an MHC structure by allele.

WHY HLA IS ITS OWN MODULE (design rationale, kept per the brief):
  Antibody recognition = a paratope reads a native antigen surface directly.
  HLA/MHC recognition = the molecule presents a short *processed peptide* inside
  a groove to a T-cell receptor (the peptide-MHC-TCR ternary system). Different
  recognition geometry, different downstream biology (transplant matching, drug
  hypersensitivity, autoimmune/disease association). Hence a distinct module,
  not a branch of the antibody layer.

What this module does (dependency-free):
  * Detects whether a loaded structure IS an MHC molecule (class I = polymorphic
    heavy chain + beta-2-microglobulin; class II = alpha + beta chains).
  * Annotates the peptide-binding groove from a FIXED reference set of
    groove-lining positions (pockets A-F for class I) — NOT blind LIGSITE, which
    would mis-find an enzyme-style cavity. Reference positions are mapped onto the
    loaded structure's numbering by sequence alignment (reusing the NW aligner).
  * Compares how THIS allele's groove differs from a class-I reference: which
    pocket-lining positions carry allele-specific residues and how that shifts
    charge / hydrophobicity / size (reusing the Pockets residue-property sets).
  * Cross-references a curated HLA <-> drug hypersensitivity table (a lookup
    table, NOT a prediction model). RESEARCH ONLY — clinically adjacent.

Peptide-into-groove docking (brief 3.5) is DEFERRED behind an "experimental"
flag (see ``PEPTIDE_DOCKING_AVAILABLE``) to stay within the synchronous budget.

Allele input assumes the user KNOWS/SPECIFIES the allele (typed or implied by the
loaded PDB). HLA typing from raw sequencing reads is explicitly out of scope.
"""

from __future__ import annotations

import re

from .evolution import AA3TO1, _nw_align
from .pockets import AROMATIC_RES, CHARGED_RES, HYDROPHOBIC_RES, POLAR_RES

PEPTIDE_DOCKING_AVAILABLE = False  # brief 3.5 — deferred/experimental

# ---------------------------------------------------------------------------
# Reference sequences (mature protein numbering).
# ---------------------------------------------------------------------------

# HLA-A*02:01 mature alpha chain — the canonical class-I reference. Groove
# positions below are in this (standard) numbering; the structure's heavy chain
# is aligned onto it so any author-numbering offset is handled.
CLASS_I_REF = (
    "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQ"
    "THRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTK"
    "HKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT"
    "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP"
)

# Human beta-2-microglobulin (mature) — highly conserved; class-I marker chain.
B2M_REF = (
    "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTP"
    "TEKDEYACRVNHVTLSQPKIVKWDRDM"
)

# HLA-DRA (alpha) and a DRB1 (beta) domain reference — for best-effort class-II
# detection and groove annotation (class II is polymorphic mainly on the beta
# chain).
CLASS_II_ALPHA_REF = (
    "IKEEHVIIQAEFYLNPDQSGEFMFDFDGDEIFHVDMAKKETVWRLEEFGRFASFEAQGALANIAVDKANLEI"
    "MTKRSNYTPITN"
)
CLASS_II_BETA_REF = (
    "GDTRPRFLWQLKFECHFFNGTERVRLLERCIYNQEESVRFDSDVGEYRAVTELGRPDAEYWNSQKDLLEQRR"
    "AAVDTYCRHNYGVGESFTVQRR"
)

# Class-I peptide-binding groove: pocket -> lining residue positions
# (HLA class-I mature numbering; Saper, Bjorkman & Wiley 1991, J Mol Biol 219:277;
# Madden 1995, Annu Rev Immunol 13:587). Pockets A and F anchor the peptide
# termini; B-E accommodate side chains along the peptide.
CLASS_I_POCKETS = {
    "A": [5, 59, 63, 66, 159, 163, 167, 171],
    "B": [7, 9, 24, 34, 45, 63, 66, 67, 70, 99],
    "C": [9, 70, 73, 74, 97],
    "D": [99, 114, 155, 156, 159, 160],
    "E": [97, 114, 133, 147, 152],
    "F": [77, 80, 81, 84, 95, 116, 123, 143, 146, 147],
}

# Class-II groove-lining positions (best-effort; beta1 dominates specificity).
CLASS_II_POCKETS = {
    "beta1": [9, 11, 13, 26, 28, 30, 37, 47, 57, 60, 61, 70, 71, 74, 78, 81, 85, 86, 89, 90],
    "alpha1": [9, 11, 22, 24, 31, 43, 52, 53, 62, 65, 66, 68, 69, 72, 76],
}

# Approximate side-chain volumes (Å³) for size-change description.
_VOLUME = {
    "GLY": 60, "ALA": 89, "SER": 89, "CYS": 109, "ASP": 111, "PRO": 113,
    "ASN": 114, "THR": 116, "GLU": 138, "VAL": 140, "GLN": 144, "HIS": 153,
    "MET": 163, "ILE": 167, "LEU": 167, "LYS": 169, "ARG": 173, "PHE": 190,
    "TYR": 194, "TRP": 228,
}
_POS_RES = {"LYS", "ARG", "HIS"}
_NEG_RES = {"ASP", "GLU"}


def _charge(res_name):
    if res_name in _POS_RES:
        return 1
    if res_name in _NEG_RES:
        return -1
    return 0


def _property_class(res_name):
    if res_name in CHARGED_RES:
        return "charged"
    if res_name in HYDROPHOBIC_RES:
        return "hydrophobic"
    if res_name in POLAR_RES:
        return "polar"
    return "other"


# ---------------------------------------------------------------------------
# Allele nomenclature
# ---------------------------------------------------------------------------

_ALLELE_RE = re.compile(
    r"^(?:HLA[-\s]?)?"
    r"(A|B|C|E|F|G|DRA|DRB[1-9]|DQA1|DQB1|DPA1|DPB1)"
    r"[\*\s]?"
    r"(\d{2,4})"
    r"(?::(\d{2,3}))?",
    re.IGNORECASE,
)


def normalize_allele(text):
    """Normalize an HLA allele string to 'LOCUS*FF:FF' (2-field), or None."""
    if not text:
        return None
    m = _ALLELE_RE.match(text.strip())
    if not m:
        return None
    locus = m.group(1).upper()
    f1 = m.group(2)
    f2 = m.group(3)
    if f2 is None:
        # No colon: a 4-digit run is the two fields concatenated (e.g. 'B5701').
        if len(f1) >= 4:
            f1, f2 = f1[:2], f1[2:4]
        else:
            return f"{locus}*{f1}"
    return f"{locus}*{f1}:{f2}"


def allele_class(allele):
    """Return 'I' or 'II' from a normalized allele's locus."""
    if not allele:
        return None
    locus = allele.split("*")[0]
    if locus in ("A", "B", "C", "E", "F", "G"):
        return "I"
    if locus.startswith(("DR", "DQ", "DP")):
        return "II"
    return None


# ---------------------------------------------------------------------------
# HLA <-> drug hypersensitivity lookup table (curated; NOT a prediction model)
# ---------------------------------------------------------------------------

DRUG_ASSOCIATIONS = {
    "B*57:01": [{
        "drug": "abacavir",
        "reaction": "abacavir hypersensitivity syndrome",
        "citation": "Mallal et al. 2008, NEJM 358:568 (PREDICT-1)",
    }],
    "B*15:02": [{
        "drug": "carbamazepine",
        "reaction": "Stevens-Johnson syndrome / toxic epidermal necrolysis",
        "citation": "Chung et al. 2004, Nature 428:486",
    }],
    "A*31:01": [{
        "drug": "carbamazepine",
        "reaction": "carbamazepine-induced hypersensitivity (incl. DRESS)",
        "citation": "McCormack et al. 2011, NEJM 364:1134",
    }],
    "B*58:01": [{
        "drug": "allopurinol",
        "reaction": "severe cutaneous adverse reaction (SJS/TEN, DRESS)",
        "citation": "Hung et al. 2005, PNAS 102:4134",
    }],
    "B*13:01": [{
        "drug": "dapsone",
        "reaction": "dapsone hypersensitivity syndrome",
        "citation": "Zhang et al. 2013, NEJM 369:1620",
    }],
    "A*32:01": [{
        "drug": "vancomycin",
        "reaction": "vancomycin-induced DRESS",
        "citation": "Konvinse et al. 2019, J Allergy Clin Immunol 144:183",
    }],
    "B*57:01 ": [],  # guard against accidental trailing space keys
}
# Drop any malformed keys.
DRUG_ASSOCIATIONS = {k.strip(): v for k, v in DRUG_ASSOCIATIONS.items() if v}


def drug_hits(allele):
    """Return curated drug-hypersensitivity associations for a normalized allele."""
    if not allele:
        return []
    return list(DRUG_ASSOCIATIONS.get(allele, []))


# ---------------------------------------------------------------------------
# Structure sequence helpers
# ---------------------------------------------------------------------------

def _chain_sequences(structure):
    by_chain = {}
    for a in structure.protein_atoms:
        d = by_chain.setdefault(a.chain, {})
        if a.res_seq not in d:
            d[a.res_seq] = a.res_name
    out = {}
    for chain, residues in by_chain.items():
        items = sorted(residues.items())
        seq = "".join(AA3TO1.get(rn, "X") for _, rn in items)
        out[chain] = (seq, items)  # items: [(res_seq, res_name)]
    return out


def _identity_to(seq, ref, window=320):
    """NW-align seq (capped) to ref; return (identity_over_ref_alignment, pairs)."""
    s = seq[:window]
    pairs = _nw_align(s, ref)
    aligned = matches = 0
    for ti, rj in pairs:
        if ti is None or rj is None:
            continue
        aligned += 1
        if s[ti] == ref[rj]:
            matches += 1
    return (matches / aligned if aligned else 0.0), pairs


# ---------------------------------------------------------------------------
# Detection (cheap; runs inline on structure load)
# ---------------------------------------------------------------------------

def detect(structure):
    """Classify a structure as MHC class I / II / not-MHC from chain sequences."""
    seqs = _chain_sequences(structure)
    heavy = b2m = alpha = beta = None
    heavy_id = b2m_id = alpha_id = beta_id = 0.0
    for chain, (seq, _items) in seqs.items():
        if len(seq) < 60:
            continue
        idc1, _ = _identity_to(seq, CLASS_I_REF)
        if idc1 > heavy_id:
            heavy_id, heavy = idc1, chain
        idb2m, _ = _identity_to(seq, B2M_REF)
        if idb2m > b2m_id:
            b2m_id, b2m = idb2m, chain
        ida, _ = _identity_to(seq, CLASS_II_ALPHA_REF)
        if ida > alpha_id:
            alpha_id, alpha = ida, chain
        idb, _ = _identity_to(seq, CLASS_II_BETA_REF)
        if idb > beta_id:
            beta_id, beta = idb, chain

    is_class_i = heavy_id >= 0.55 and heavy is not None
    is_b2m = b2m_id >= 0.6 and b2m is not None and b2m != heavy
    # Class II: distinct alpha and beta chains both matching, and NOT a stronger
    # class-I heavy match.
    is_class_ii = (
        alpha_id >= 0.55 and beta_id >= 0.55 and alpha != beta
        and not is_class_i
    )

    if is_class_i:
        return {
            "is_mhc": True, "mhc_class": "I",
            "heavy_chain": heavy, "b2m_chain": b2m if is_b2m else None,
            "heavy_identity": round(heavy_id, 3),
            "b2m_identity": round(b2m_id, 3) if is_b2m else None,
        }
    if is_class_ii:
        return {
            "is_mhc": True, "mhc_class": "II",
            "alpha_chain": alpha, "beta_chain": beta,
            "alpha_identity": round(alpha_id, 3), "beta_identity": round(beta_id, 3),
        }
    return {"is_mhc": False}


# ---------------------------------------------------------------------------
# Groove annotation + allele-specific comparison
# ---------------------------------------------------------------------------

def _map_ref_positions(seq, items, ref):
    """Map reference (1-based) positions -> {pos: (res_seq, res_name, struct_aa, ref_aa)}."""
    identity, pairs = _identity_to(seq, ref, window=len(seq))
    out = {}
    for ti, rj in pairs:
        if ti is None or rj is None:
            continue
        res_seq, res_name = items[ti]
        out[rj + 1] = (res_seq, res_name, seq[ti], ref[rj])
    return out, identity


def _annotate_groove(seq, items, ref, pockets_def, ref_seq):
    """Build groove residue annotations + allele-vs-reference differences."""
    pos_map, identity = _map_ref_positions(seq, items, ref)
    residues = []
    differences = []
    pocket_summ = {}
    seen = {}
    for pocket, positions in pockets_def.items():
        summ = {"pocket": pocket, "positions": len(positions), "mapped": 0,
                "polymorphic": 0, "charge_shift": 0, "detail": []}
        for pos in positions:
            entry = pos_map.get(pos)
            if entry is None:
                continue
            summ["mapped"] += 1
            res_seq, res_name, struct_aa, _ = entry
            ref_res = _three_from_one(ref_seq[pos - 1]) if pos - 1 < len(ref_seq) else None
            key = (res_seq,)
            rec = seen.get(key)
            if rec is None:
                rec = {"res_seq": res_seq, "res_name": res_name,
                       "pockets": [], "ref_res": ref_res,
                       "polymorphic": res_name != ref_res if ref_res else False}
                seen[key] = rec
                residues.append(rec)
            if pocket not in rec["pockets"]:
                rec["pockets"].append(pocket)

            if ref_res and res_name != ref_res:
                summ["polymorphic"] += 1
                dq = _charge(res_name) - _charge(ref_res)
                if dq:
                    summ["charge_shift"] += dq
                diff = {
                    "res_seq": res_seq, "pocket": pocket,
                    "reference": ref_res, "allele": res_name,
                    "change": _describe_change(ref_res, res_name),
                }
                differences.append(diff)
                summ["detail"].append(f"{ref_res}{pos}{res_name}")
        pocket_summ[pocket] = summ

    # De-duplicate differences that recur across shared positions.
    uniq = {}
    for d in differences:
        uniq.setdefault((d["res_seq"], d["allele"]), d)
    return {
        "reference_identity": round(identity, 3),
        "groove_residues": sorted(residues, key=lambda r: r["res_seq"]),
        "differences": list(uniq.values()),
        "pocket_summary": list(pocket_summ.values()),
    }


# Canonical 1->3 letter map (standard 20; AA3TO1 aliases MSE->M etc. so invert
# explicitly to stay unambiguous).
_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def _three_from_one(aa):
    return _ONE_TO_THREE.get(aa)


def _describe_change(ref_res, allele_res):
    tags = []
    dq = _charge(allele_res) - _charge(ref_res)
    if dq > 0:
        tags.append("more positive")
    elif dq < 0:
        tags.append("more negative")
    pc_ref, pc_alt = _property_class(ref_res), _property_class(allele_res)
    if pc_ref != pc_alt:
        tags.append(f"{pc_ref}→{pc_alt}")
    dv = _VOLUME.get(allele_res, 0) - _VOLUME.get(ref_res, 0)
    if abs(dv) >= 40:
        tags.append("larger" if dv > 0 else "smaller")
    if allele_res in AROMATIC_RES and ref_res not in AROMATIC_RES:
        tags.append("gains aromatic")
    return ", ".join(tags) or "conservative"


def analyze(structure, allele=None):
    """Full HLA analysis: detection, groove annotation, allele diffs, drug hits."""
    det = detect(structure)
    norm_allele = normalize_allele(allele) if allele else None
    cls = det.get("mhc_class") or allele_class(norm_allele)

    result = {
        "is_mhc": det.get("is_mhc", False),
        "detection": det,
        "allele": norm_allele,
        "allele_input": allele,
        "mhc_class": cls,
        "peptide_docking_available": PEPTIDE_DOCKING_AVAILABLE,
    }

    seqs = _chain_sequences(structure)

    if cls == "I" and det.get("heavy_chain") in seqs:
        seq, items = seqs[det["heavy_chain"]]
        groove = _annotate_groove(seq, items, CLASS_I_REF, CLASS_I_POCKETS, CLASS_I_REF)
        groove["chain"] = det["heavy_chain"]
        groove["reference_allele"] = "A*02:01"
        result["groove"] = groove
    elif cls == "II" and det.get("beta_chain") in seqs:
        # Best-effort: annotate the beta chain groove (specificity-dominant).
        seq, items = seqs[det["beta_chain"]]
        groove = _annotate_groove(
            seq, items, CLASS_II_BETA_REF,
            {"beta1": CLASS_II_POCKETS["beta1"]}, CLASS_II_BETA_REF,
        )
        groove["chain"] = det["beta_chain"]
        groove["reference_allele"] = "DRB1 reference"
        groove["approximate"] = True
        result["groove"] = groove

    if norm_allele:
        result["drug_associations"] = drug_hits(norm_allele)

    return result
