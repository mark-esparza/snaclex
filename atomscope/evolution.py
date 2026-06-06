"""Evolutionary conservation analysis from a protein-family alignment.

Pipeline (dependency-free):
  1. UniProt accession -> Pfam family (InterPro API) -> family alignment.
  2. Per-column Shannon-entropy conservation + consensus residue.
  3. Extract the loaded structure's sequence (with residue numbering).
  4. Needleman-Wunsch align the target to the family consensus to map each
     structure residue onto an alignment column, giving per-residue conservation.
  5. Score each detected pocket by the mean conservation of its lining residues.

Conservation here is a family-MSA signal (how invariant a position is across
homologs), not a phylogenetic reconstruction. Research-only.
"""

from __future__ import annotations

import gzip
import math
import urllib.error
import urllib.parse
import urllib.request

from .http_util import FetchError, fetch_bytes, fetch_json

MAX_ALIGN_SEQS = 1200  # cap on homologs used (deep enough for coevolution, bounded)
COEVO_CONF_MIN = 0.6   # min fraction of top pairs that are spatial contacts to trust the signal
COEVO_CONTACT_A = 12.0  # CA-CA distance counted as a structural contact
_UA = "AtomScope/0.1 (research tool; +local)"

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}
_AA = "ACDEFGHIKLMNPQRSTVWY"
_LOG20 = math.log2(20)


# ---------------- alignment acquisition ----------------

def pfam_for_uniprot(accession: str) -> dict | None:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}/"
    try:
        data = fetch_json(url)
    except FetchError:
        return None
    results = data.get("results") or []
    if not results:
        return None
    meta = results[0].get("metadata", {})
    return {"pfam": meta.get("accession"), "name": meta.get("name")}


def fetch_family_alignment(pfam_id: str) -> list[tuple]:
    """Return [(name, aligned_seq)] homologs for a Pfam family.

    Streams the *full* alignment (deep enough for coevolution) but stops after
    MAX_ALIGN_SEQS sequences to bound time/memory. Falls back to the small seed
    alignment if streaming the full one fails.
    """
    try:
        aln = _stream_alignment(pfam_id, "alignment:full", MAX_ALIGN_SEQS)
        if len(aln) >= 50:
            return _modal_width(aln)
    except (FetchError, OSError, urllib.error.URLError, TimeoutError):
        pass
    try:
        raw = fetch_bytes(
            f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/"
            "?annotation=alignment:seed"
        )
        text = gzip.decompress(raw).decode("utf-8", "replace")
        return _modal_width(_parse_stockholm(text.splitlines()))
    except (FetchError, OSError, EOFError):
        return []


def _stream_alignment(pfam_id: str, kind: str, max_seqs: int) -> list[tuple]:
    url = (
        f"https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam_id}/"
        f"?annotation={urllib.parse.quote(kind)}"
    )
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    resp = urllib.request.urlopen(req, timeout=60)
    rows = []
    try:
        gz = gzip.GzipFile(fileobj=resp)
        for raw_line in gz:
            line = raw_line.decode("utf-8", "replace").rstrip("\n")
            if not line or line[0] == "#" or line.startswith("//"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            rows.append((parts[0], parts[1].strip().upper().replace(".", "-")))
            if len(rows) >= max_seqs:
                break
    finally:
        resp.close()
    return rows


def _parse_stockholm(lines) -> list[tuple]:
    rows: dict[str, list] = {}
    order: list[str] = []
    for line in lines:
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        name, seq = parts
        if name not in rows:
            rows[name] = []
            order.append(name)
        rows[name].append(seq.strip())
    return [(n, "".join(rows[n]).upper().replace(".", "-")) for n in order]


def _modal_width(aln: list[tuple]) -> list[tuple]:
    if not aln:
        return []
    from collections import Counter

    width = Counter(len(s) for _, s in aln).most_common(1)[0][0]
    return [(n, s) for (n, s) in aln if len(s) == width]


# ---------------- conservation ----------------

def column_conservation(alignment: list[tuple]):
    """Return (conservation[], consensus[], occupancy[]) per alignment column."""
    if not alignment:
        return [], [], []
    width = len(alignment[0][1])
    n = len(alignment)
    conservation, consensus, occupancy = [], [], []
    for col in range(width):
        counts: dict[str, int] = {}
        nongap = 0
        for _, seq in alignment:
            c = seq[col]
            if c == "-" or c not in _AA:
                continue
            counts[c] = counts.get(c, 0) + 1
            nongap += 1
        if nongap == 0:
            conservation.append(0.0)
            consensus.append("-")
            occupancy.append(0.0)
            continue
        entropy = 0.0
        for c, k in counts.items():
            p = k / nongap
            entropy -= p * math.log2(p)
        cons = 1.0 - entropy / _LOG20
        conservation.append(max(0.0, min(cons, 1.0)))
        consensus.append(max(counts, key=counts.get))
        occupancy.append(nongap / n)
    return conservation, consensus, occupancy


# ---------------- target sequence from structure ----------------

def structure_sequence(structure):
    """Pick the longest protein chain; return (one_letter_seq, [(chain,res_seq,res_name)])."""
    by_chain: dict[str, dict] = {}
    for a in structure.protein_atoms:
        by_chain.setdefault(a.chain, {})
        if a.res_seq not in by_chain[a.chain]:
            by_chain[a.chain][a.res_seq] = a.res_name
    if not by_chain:
        return "", []
    chain = max(by_chain, key=lambda c: len(by_chain[c]))
    residues = sorted(by_chain[chain].items())  # (res_seq, res_name)
    seq = []
    keys = []
    for res_seq, res_name in residues:
        seq.append(AA3TO1.get(res_name, "X"))
        keys.append((chain, res_seq, res_name))
    return "".join(seq), keys


# ---------------- Needleman-Wunsch (global) ----------------

def _nw_align(a: str, b: str):
    """Global alignment; returns list of (i|None, j|None) index pairs."""
    GAP = -4
    n, m = len(a), len(b)
    # Score matrix rows kept as lists; traceback via pointer matrix.
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = i * GAP
    for j in range(1, m + 1):
        score[0][j] = j * GAP
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = score[i], score[i - 1]
        for j in range(1, m + 1):
            sub = 5 if ai == b[j - 1] else -1
            row[j] = max(prev[j - 1] + sub, prev[j] + GAP, row[j - 1] + GAP)
    # Traceback.
    pairs = []
    i, j = n, m
    while i > 0 and j > 0:
        cur = score[i][j]
        sub = 5 if a[i - 1] == b[j - 1] else -1
        if cur == score[i - 1][j - 1] + sub:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif cur == score[i - 1][j] + GAP:
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    while i > 0:
        pairs.append((i - 1, None)); i -= 1
    while j > 0:
        pairs.append((None, j - 1)); j -= 1
    pairs.reverse()
    return pairs


# ---------------- orchestration ----------------

def analyze(structure, uniprots: list[str]) -> dict | None:
    """Full conservation analysis for a loaded structure. None if unavailable."""
    fam = None
    acc_used = None
    for acc in uniprots or []:
        fam = pfam_for_uniprot(acc)
        if fam and fam.get("pfam"):
            acc_used = acc
            break
    if not fam or not fam.get("pfam"):
        return None

    try:
        alignment = fetch_family_alignment(fam["pfam"])
    except FetchError:
        return None
    if len(alignment) < 5:
        return None

    cons, consensus, occ = column_conservation(alignment)
    # Core columns = reasonably occupied positions.
    core_cols = [i for i in range(len(cons)) if occ[i] >= 0.5]
    core_consensus = "".join(consensus[i] for i in core_cols)
    core_cons = [cons[i] for i in core_cols]
    if not core_consensus:
        return None

    target_seq, keys = structure_sequence(structure)
    if not target_seq:
        return None

    pairs = _nw_align(target_seq, core_consensus)
    residues = []
    mapped = 0
    cons_by_key: dict[tuple, float] = {}
    col_to_residue: dict[int, tuple] = {}
    for ti, cj in pairs:
        if ti is None:
            continue
        chain, res_seq, res_name = keys[ti]
        if cj is None:
            residues.append(
                {"chain": chain, "res_seq": res_seq, "res_name": res_name,
                 "conservation": None}
            )
            continue
        c = round(core_cons[cj], 3)
        mapped += 1
        cons_by_key[(chain, res_seq)] = c
        col_to_residue[cj] = (chain, res_seq, res_name)
        target_aa = AA3TO1.get(res_name, "X")
        cons_aa = core_consensus[cj]
        residues.append(
            {"chain": chain, "res_seq": res_seq, "res_name": res_name,
             "conservation": c, "consensus": cons_aa, "target_aa": target_aa,
             "divergent": target_aa != cons_aa}
        )

    ranked = sorted(
        (r for r in residues if r["conservation"] is not None),
        key=lambda r: -r["conservation"],
    )
    top_conserved = [
        {"res": f"{r['res_name']}{r['res_seq']}", "chain": r["chain"],
         "conservation": r["conservation"]}
        for r in ranked[:15]
    ]

    # Ancestral-divergence (ASR-lite): where this protein departs from the
    # conserved family consensus. A change at an otherwise-conserved column is a
    # lineage-specific "derived" substitution — often specificity-determining.
    mapped_cons = [r for r in residues if r["conservation"] is not None]
    matches = sum(1 for r in mapped_cons if not r["divergent"])
    consensus_identity = round(matches / len(mapped_cons), 3) if mapped_cons else 0.0
    notable = [
        r for r in mapped_cons if r["divergent"] and r["conservation"] >= 0.4
    ]
    notable.sort(key=lambda r: -r["conservation"])
    divergent_residues = [
        {"res": f"{r['res_name']}{r['res_seq']}", "chain": r["chain"],
         "from": r["consensus"], "to": r["target_aa"],
         "conservation": r["conservation"]}
        for r in notable[:20]
    ]
    divergent_keys = {(r["chain"], r["res_seq"]) for r in notable}

    # Evolutionary coupling (APC-corrected mutual information) -> co-evolving
    # residue pairs and coupling hubs, mapped back onto the structure.
    ca = _ca_coords(structure)
    coupling_pairs, coupling_hubs, hub_keys, coupling_conf = _coevolution(
        alignment, core_cols, col_to_residue, cons_by_key, ca
    )

    return {
        "pfam": fam["pfam"],
        "family_name": fam.get("name"),
        "uniprot": acc_used,
        "n_sequences": len(alignment),
        "target_length": len(target_seq),
        "mapped_residues": mapped,
        "coverage": round(mapped / len(target_seq), 2) if target_seq else 0,
        "residues": residues,
        "top_conserved": top_conserved,
        "coupling_pairs": coupling_pairs,
        "coupling_hubs": coupling_hubs,
        "coupling_confidence": coupling_conf,
        "coupling_reliable": coupling_conf >= COEVO_CONF_MIN and bool(coupling_pairs),
        "consensus_identity": consensus_identity,
        "divergent_residues": divergent_residues,
        "_cons_by_key": cons_by_key,      # internal: for pocket scoring
        "_hub_keys": hub_keys,            # internal: for pocket allosteric flag
        "_divergent_keys": divergent_keys,  # internal: for pocket divergence
    }


def _ca_coords(structure) -> dict:
    """Representative coordinate per residue (CA preferred)."""
    ca: dict[tuple, tuple] = {}
    for a in structure.protein_atoms:
        key = (a.chain, a.res_seq)
        if a.name == "CA":
            ca[key] = (a.x, a.y, a.z)
        elif key not in ca:
            ca[key] = (a.x, a.y, a.z)
    return ca


def _coevolution(alignment, core_cols, col_to_residue, cons_by_key, ca,
                 top_n=40, min_sep=5):
    """APC-corrected mutual information between alignment columns.

    Returns (pairs, hubs, hub_key_set). Only columns that map to a structure
    residue are considered, so every reported pair is structurally locatable.
    """
    seqs = [s for _, s in alignment]
    n = len(seqs)
    if n < 50:  # too shallow for a trustworthy coupling signal
        return [], [], set(), 0.0
    sym = {a: i for i, a in enumerate(_AA)}

    # Encode mapped columns; keep only well-occupied, *variable* columns —
    # invariant/conserved columns cannot co-vary, and skipping them is the main
    # speed-up (and de-noises the signal).
    enc = {}
    keep = []
    half = 0.5 * n
    log2 = math.log2
    for ci in range(len(core_cols)):
        if ci not in col_to_residue:
            continue
        oc = core_cols[ci]
        col = [sym.get(s[oc], 20) for s in seqs]
        cnt = {}
        m = 0
        for v in col:
            if v != 20:
                cnt[v] = cnt.get(v, 0) + 1
                m += 1
        if m < half:
            continue
        entropy = 0.0
        for c in cnt.values():
            p = c / m
            entropy -= p * log2(p)
        if entropy < 0.5:
            continue
        enc[ci] = col
        keep.append(ci)

    L = len(keep)
    if L < 4:
        return [], [], set(), 0.0

    mi = {}
    mi_row_sum = {c: 0.0 for c in keep}
    for ii in range(L):
        ci = keep[ii]
        col_i = enc[ci]
        oci = core_cols[ci]
        for jj in range(ii + 1, L):
            cj = keep[jj]
            if abs(core_cols[cj] - oci) < min_sep:
                continue
            col_j = enc[cj]
            joint = [0] * 441
            fa = [0] * 21
            fb = [0] * 21
            m = 0
            for k in range(n):
                a = col_i[k]
                if a == 20:
                    continue
                b = col_j[k]
                if b == 20:
                    continue
                joint[a * 21 + b] += 1
                fa[a] += 1
                fb[b] += 1
                m += 1
            if m < half:
                continue
            present_a = [a for a in range(21) if fa[a]]
            present_b = [b for b in range(21) if fb[b]]
            val = 0.0
            for a in present_a:
                pa = fa[a] / m
                base = a * 21
                for b in present_b:
                    c = joint[base + b]
                    if c:
                        pab = c / m
                        val += pab * log2(pab / (pa * (fb[b] / m)))
            if val > 0:
                mi[(ci, cj)] = val
                mi_row_sum[ci] += val
                mi_row_sum[cj] += val

    if not mi:
        return [], [], set(), 0.0
    overall = sum(mi.values()) * 2 / (L * (L - 1)) if L > 1 else 0.0
    col_mean = {c: mi_row_sum[c] / (L - 1) for c in keep} if L > 1 else {}

    scored = []
    for (ci, cj), val in mi.items():
        apc = (col_mean[ci] * col_mean[cj] / overall) if overall > 0 else 0.0
        scored.append((ci, cj, val - apc))
    scored.sort(key=lambda x: -x[2])

    pairs = []
    degree: dict[tuple, float] = {}
    for ci, cj, mip in scored[:top_n]:
        ri, rj = col_to_residue[ci], col_to_residue[cj]
        ki, kj = (ri[0], ri[1]), (rj[0], rj[1])
        xi, xj = ca.get(ki), ca.get(kj)
        dist = None
        if xi and xj:
            dist = round(
                math.dist(xi, xj), 1
            )
        pairs.append({
            "res_i": f"{ri[2]}{ri[1]}", "chain_i": ri[0],
            "res_j": f"{rj[2]}{rj[1]}", "chain_j": rj[0],
            "mip": round(mip, 3),
            "distance_A": dist,
            "xyz_i": [round(v, 3) for v in xi] if xi else None,
            "xyz_j": [round(v, 3) for v in xj] if xj else None,
        })
        degree[ki] = degree.get(ki, 0) + 1
        degree[kj] = degree.get(kj, 0) + 1

    # Self-validation: what fraction of top pairs are actual spatial contacts?
    # This is family-dependent; a shallow/divergent MSA yields a noisy signal we
    # should NOT present as real coupling.
    dists = [p["distance_A"] for p in pairs if p["distance_A"] is not None]
    confidence = (
        sum(1 for d in dists if d <= COEVO_CONTACT_A) / len(dists) if dists else 0.0
    )
    if confidence < COEVO_CONF_MIN:
        return [], [], set(), round(confidence, 2)

    name_of = {(r[0], r[1]): r[2] for r in col_to_residue.values()}
    hubs = sorted(degree.items(), key=lambda kv: -kv[1])
    coupling_hubs = []
    for (chain, res_seq), deg in hubs[:12]:
        coupling_hubs.append({
            "res": f"{name_of.get((chain, res_seq), '')}{res_seq}", "chain": chain,
            "degree": deg, "conservation": cons_by_key.get((chain, res_seq)),
        })
    # A "hub" for the allosteric flag is a genuine network node (degree >= 3),
    # so the flag stays rare and meaningful (not every coupled residue).
    hub_keys = {k for k, deg in degree.items() if deg >= 3}
    return pairs, coupling_hubs, hub_keys, round(confidence, 2)


def annotate_pockets(evo: dict, pockets: list) -> list:
    """Add mean conservation + coupling-hub count + an evolutionary label.

    Conserved pocket -> likely functional; coupled-but-not-conserved pocket ->
    candidate allosteric/cryptic site (the co-evolution signal the literature
    associates with hidden control pockets).
    """
    cons_by_key = evo.get("_cons_by_key", {})
    hub_keys = evo.get("_hub_keys", set())
    divergent_keys = evo.get("_divergent_keys", set())
    out = []
    for p in pockets:
        vals = []
        hubs = 0
        divergent = 0
        for r in p.get("lining_residues", []):
            key = (r["chain"], r["res_seq"])
            c = cons_by_key.get(key)
            if c is not None:
                vals.append(c)
            if key in hub_keys:
                hubs += 1
            if key in divergent_keys:
                divergent += 1
        if vals:
            mean = sum(vals) / len(vals)
            enriched = hubs >= 3 and len(vals) and hubs >= 0.25 * len(vals)
            if enriched and mean < 0.42:
                label = "candidate allosteric (coupling-enriched, not conserved)"
            elif mean >= 0.55:
                label = "evolutionarily conserved (likely functional)"
            elif mean >= 0.35:
                label = "moderately conserved"
            else:
                label = "variable"
        else:
            mean = None
            label = "unmapped"
        out.append(
            {
                "index": p["index"],
                "tier": p.get("tier"),
                "volume_A3": p.get("volume_A3"),
                "mean_conservation": round(mean, 3) if mean is not None else None,
                "conserved_residues": len(vals),
                "coupling_hubs": hubs,
                "divergent_lining": divergent,
                "specificity_candidate": divergent >= 2 and len(vals) and divergent >= 0.2 * len(vals),
                "label": label,
            }
        )
    return out
