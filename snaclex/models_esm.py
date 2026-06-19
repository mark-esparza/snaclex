"""EvolutionaryScale ESM model integration (optional, env-gated).

Two capabilities, both OFF unless ``ESM_API_KEY`` is set:

  * ``fold_sequence(seq)`` — a predicted structure (PDB) from a bare sequence,
    so SnaCleX can analyze an HLA allele's protein (or any protein) with no
    experimental structure. The result flows through the normal structure
    pipeline (pockets / interactions / groove).
  * ``score_variants(seq, variants)`` — a sequence-model log-likelihood-ratio
    per substitution: a transparent variant-effect "predictability" signal that
    sits ALONGSIDE the structural overlay, never replacing it.

Both call EvolutionaryScale's Forge API (``forge.evolutionaryscale.ai``). The
network call is isolated behind an injectable seam (``scorer`` / ``folder``), so
orchestration is unit-testable offline and the build runs fully without a token
— RCSB + structural analysis remain the default path. Outputs are predictions,
labelled as such, and are never clinical truth.

Note: the concrete Forge request contract is deployment-specific; the request
URLs are taken from ``SNACLEX_ESM_SCORE_URL`` / ``SNACLEX_ESM_FOLD_URL`` so the
integration can be pointed at the official SDK-backed endpoint (or a proxy)
without code changes. When unset, the model features report unavailable.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

from .variants import parse_variant

_DEFAULT_BASE = "https://forge.evolutionaryscale.ai"
FOLD_MODEL = "esm3-open"
SCORE_MODEL = "esmc-600m"
MAX_FOLD_RESIDUES = 1024

_AA = set("ACDEFGHIKLMNPQRSTVWY")


def token() -> str:
    return os.environ.get("ESM_API_KEY") or ""


def available() -> bool:
    """True only when a Forge token is configured."""
    return bool(token())


def config() -> dict:
    return {
        "available": available(),
        "fold_model": FOLD_MODEL,
        "score_model": SCORE_MODEL,
        "base_url": os.environ.get("SNACLEX_ESM_API") or _DEFAULT_BASE,
        "max_fold_residues": MAX_FOLD_RESIDUES,
        "provider": "EvolutionaryScale Forge",
    }


def _clean_seq(seq: str) -> str:
    return "".join(c for c in (seq or "").upper() if c in _AA)


def _post_json(url: str, payload: dict, timeout: int = 120) -> dict:
    headers = {
        "Content-Type": "application/json",
        "User-Agent": "SnaCleX/0.1 (research tool)",
    }
    if token():
        headers["Authorization"] = f"Bearer {token()}"
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(url, data=data, headers=headers, method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8", errors="replace"))


# ---------------- variant-effect scoring ----------------

def _effect_label(llr: float) -> str:
    """Bucket a log-likelihood ratio (mut vs wt) into a coarse effect label."""
    if llr <= -6:
        return "strongly disfavored"
    if llr <= -2:
        return "disfavored"
    if llr < 2:
        return "tolerated"
    return "favored"


def _forge_score(seq: str, pos: int, wt: str, mut: str) -> float:
    """Call Forge for a masked-marginal LLR. Raises if not configured/reachable."""
    url = os.environ.get("SNACLEX_ESM_SCORE_URL")
    if not url:
        raise RuntimeError("Forge scoring endpoint not configured (SNACLEX_ESM_SCORE_URL)")
    data = _post_json(url, {"model": SCORE_MODEL, "sequence": seq,
                            "position": pos, "wt": wt, "mut": mut})
    return float(data["llr"])


def score_variants(seq: str, variants, scorer=None) -> dict:
    """Score substitutions by sequence-model log-likelihood ratio (offline-safe).

    ``scorer(seq, pos, wt, mut) -> float`` is injectable for testing; the default
    calls Forge. Returns ``available: False`` (without raising) when no token /
    scorer is present.
    """
    if scorer is None and not available():
        return {"available": False, "reason": "ESM_API_KEY not set"}
    seq = _clean_seq(seq)
    scorer = scorer or _forge_score

    results = []
    for v in variants:
        parsed = parse_variant(v) if isinstance(v, str) else v
        if not parsed or not parsed.get("wt") or parsed.get("mut") in (None, "*", "="):
            results.append({"input": v if isinstance(v, str) else parsed,
                            "scored": False, "reason": "not a missense substitution"})
            continue
        pos, wt, mut = parsed["position"], parsed["wt"], parsed["mut"]
        entry = {"input": parsed.get("raw"), "position": pos, "wt": wt, "mut": mut}
        if pos < 1 or pos > len(seq):
            entry.update(scored=False, reason="position outside the supplied sequence")
        elif seq[pos - 1] != wt:
            entry.update(scored=False,
                         reason=f"WT mismatch (sequence has {seq[pos - 1]} at {pos})")
        else:
            try:
                llr = float(scorer(seq, pos, wt, mut))
            except Exception as exc:  # noqa: BLE001 - degrade gracefully
                entry.update(scored=False, reason=f"scoring unavailable ({exc})")
            else:
                entry.update(scored=True, llr=round(llr, 4),
                             effect=_effect_label(llr), model=SCORE_MODEL)
        results.append(entry)

    return {
        "available": True,
        "model": SCORE_MODEL,
        "scored": results,
        "note": "Sequence-model likelihood proxy (a predictability signal), not "
                "a binding-affinity or pathogenicity prediction.",
    }


# ---------------- structure-from-sequence ----------------

def _forge_fold(seq: str) -> str:
    url = os.environ.get("SNACLEX_ESM_FOLD_URL")
    if not url:
        raise RuntimeError("Forge fold endpoint not configured (SNACLEX_ESM_FOLD_URL)")
    data = _post_json(url, {"model": FOLD_MODEL, "sequence": seq})
    pdb = data.get("pdb") or data.get("pdb_string")
    if not pdb:
        raise RuntimeError("Forge fold response contained no PDB")
    return pdb


def fold_sequence(seq: str, folder=None) -> dict:
    """Predict a structure (PDB text) from a sequence (offline-safe).

    ``folder(seq) -> pdb_text`` is injectable for testing; the default calls
    Forge. Returns ``available: False`` (without raising) on any failure.
    """
    seq = _clean_seq(seq)
    if not seq:
        return {"available": False, "reason": "empty/invalid sequence"}
    if len(seq) > MAX_FOLD_RESIDUES:
        return {"available": False,
                "reason": f"sequence too long ({len(seq)} > {MAX_FOLD_RESIDUES} residues)"}
    if folder is None and not available():
        return {"available": False, "reason": "ESM_API_KEY not set"}
    folder = folder or _forge_fold
    try:
        pdb_text = folder(seq)
    except Exception as exc:  # noqa: BLE001 - degrade gracefully
        return {"available": False, "reason": f"fold unavailable ({exc})"}
    return {
        "available": True,
        "model": FOLD_MODEL,
        "source": "esm3-predicted",
        "pdb": pdb_text,
        "length": len(seq),
        "note": "Predicted structure (ESM3); per-region confidence varies and it "
                "is not an experimental structure.",
    }
