"""Map protein variants onto a loaded structure's residues.

Given variants expressed at protein (UniProt) positions — e.g. ``R273H``,
``p.Arg273His``, or a bare position ``273`` — this aligns the UniProt canonical
sequence to the structure's own sequence (Needleman-Wunsch, reusing
``evolution._nw_align``) so each variant lands on a concrete ``(chain, residue)``
in the 3D model. The server then layers on a structural read-out (pocket /
interface / conservation) and, for genomic input, population allele frequency.

The protein-position path needs no new dependency. Genomic input
(``chr:pos ref>alt``) is gated behind Ensembl VEP and handled by the server.

Research-only: a structural localization of a variant, not a pathogenicity call.
"""

from __future__ import annotations

import re

from . import evolution
from .http_util import FetchError, fetch_text

_AA1 = set("ACDEFGHIKLMNPQRSTVWY")
# 3-letter -> 1-letter, derived from the shared evolution table.
_THREE_TO_ONE = {k.upper(): v for k, v in evolution.AA3TO1.items()}

_RE_POS = re.compile(r"^(\d+)$")
_RE_FULL = re.compile(r"^(?:p\.)?([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3}|\*|=)$")


def _norm_aa(tok: str) -> str | None:
    """Normalize a 1- or 3-letter residue token to a single letter (or * / =)."""
    t = tok.upper()
    if t in ("*", "="):
        return t
    if len(t) == 1:
        return t  # accept any single letter (X, B, Z handled downstream)
    return _THREE_TO_ONE.get(t)


def parse_variant(text: str) -> dict | None:
    """Parse a protein variant string into {raw, wt, position, mut} or None."""
    s = (text or "").strip()
    if not s:
        return None
    m = _RE_POS.match(s)
    if m:
        return {"raw": s, "wt": None, "position": int(m.group(1)), "mut": None}
    m = _RE_FULL.match(s)
    if not m:
        return None
    wt = _norm_aa(m.group(1))
    mut = _norm_aa(m.group(3))
    pos = int(m.group(2))
    if pos <= 0 or wt is None or mut is None:
        return None
    return {"raw": s, "wt": wt, "position": pos, "mut": mut}


def _parse_fasta(text: str) -> str:
    """Return the concatenated sequence from FASTA text (first record)."""
    seq = []
    for line in text.splitlines():
        if not line or line.startswith(">"):
            if seq and line.startswith(">"):
                break  # stop at the second record
            continue
        seq.append(line.strip())
    return "".join(seq).upper()


def fetch_uniprot_sequence(accession: str) -> str:
    """Fetch a UniProt canonical sequence (FASTA) and return the bare sequence."""
    acc = (accession or "").strip()
    if not acc:
        return ""
    text = fetch_text(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
    return _parse_fasta(text)


def map_positions(uniprot_seq: str, structure):
    """Map 1-based UniProt positions onto structure residues via global alignment.

    Returns ``(upos_to_key, target_seq, keys)`` where ``upos_to_key[p]`` is the
    ``(chain, res_seq, res_name)`` tuple for UniProt position ``p``.
    """
    target_seq, keys = evolution.structure_sequence(structure)
    if not target_seq or not uniprot_seq:
        return {}, target_seq, keys
    pairs = evolution._nw_align(target_seq, uniprot_seq)
    upos_to_key: dict[int, tuple] = {}
    for ti, uj in pairs:
        if ti is None or uj is None:
            continue
        upos_to_key[uj + 1] = keys[ti]
    return upos_to_key, target_seq, keys


def annotate(structure, variants, uniprot_seq: str) -> dict:
    """Locate each protein variant on the structure (no network here).

    ``variants`` is a list of strings (or pre-parsed dicts). ``uniprot_seq`` is
    the canonical sequence the variant positions refer to.
    """
    upos_to_key, target_seq, _keys = map_positions(uniprot_seq, structure)
    results = []
    for v in variants:
        raw = v if isinstance(v, str) else v.get("raw", "")
        parsed = parse_variant(raw) if isinstance(v, str) else v
        entry: dict = {"input": raw}
        if not parsed:
            entry.update(mapped=False, reason="could not parse variant")
            results.append(entry)
            continue
        entry.update(wt=parsed.get("wt"), position=parsed["position"], mut=parsed.get("mut"))
        key = upos_to_key.get(parsed["position"])
        if not key:
            entry.update(
                mapped=False,
                reason="position not covered by the structure's modeled residues",
            )
            results.append(entry)
            continue
        chain, res_seq, res_name = key
        struct_aa = evolution.AA3TO1.get(res_name, "X")
        wt = parsed.get("wt")
        entry.update(
            mapped=True,
            chain=chain,
            res_seq=res_seq,
            res_name=res_name,
            structure_aa=struct_aa,
            res_id=f"{chain}/{res_name}{res_seq}",
            wt_matches_structure=(wt is None or wt == struct_aa),
        )
        results.append(entry)

    return {
        "variants": results,
        "mapped_count": sum(1 for r in results if r.get("mapped")),
        "input_count": len(results),
        "structure_length": len(target_seq),
        "uniprot_length": len(uniprot_seq),
    }
