"""ChEMBL cross-reference: drug status, mechanism of action, and measured
bioactivity, plus a best-effort match against the loaded protein target.

Implements the docking literature's recommendation to validate predictions
against measured actives / orthogonal evidence: when a docked or inspected
chemical has curated ChEMBL pharmacology, we surface its mechanism of action,
known targets, best measured potency (Ki/IC50/Kd), and whether any of those
targets corresponds to the protein currently loaded in AtomScope.

All calls degrade gracefully to None/partial — ChEMBL coverage is uneven and
this is optional context, never a hard dependency.
"""

from __future__ import annotations

import re
import urllib.parse

from .http_util import FetchError, fetch_json

_BASE = "https://www.ebi.ac.uk/chembl/api/data"

_PHASE = {
    0: "preclinical / research compound",
    1: "Phase 1 clinical",
    2: "Phase 2 clinical",
    3: "Phase 3 clinical",
    4: "approved drug",
}

# Small in-process caches (target records and molecule resolutions are reused
# across the chemical panel, docking, and screening).
_TARGET_CACHE: dict[str, dict] = {}
_MOL_CACHE: dict[str, dict] = {}

# Tokens too generic to imply two target names are the same protein.
_STOP = {
    "human", "type", "protein", "1", "2", "3", "isoform", "subunit", "and",
    "the", "of", "alpha", "beta", "gamma", "receptor", "virus", "factor",
}


def _phase_num(value):
    try:
        return int(float(value)) if value is not None else None
    except (TypeError, ValueError):
        return None


def lookup_molecule(name: str) -> dict | None:
    """Best-effort ChEMBL lookup by name (parent-resolved). None if unavailable."""
    mol = _resolve_molecule(name)
    if not mol:
        return None
    return {
        "chembl_id": mol["chembl_id"],
        "pref_name": mol.get("pref_name"),
        "max_phase": mol.get("max_phase"),
        "development_status": _PHASE.get(mol.get("max_phase"), "unknown"),
        "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol['chembl_id']}/",
    }


def _resolve_molecule(name: str) -> dict | None:
    key = (name or "").strip().lower()
    if not key:
        return None
    if key in _MOL_CACHE:
        return _MOL_CACHE[key]
    try:
        url = f"{_BASE}/molecule/search?q={urllib.parse.quote(name)}&format=json&limit=5"
        molecules = fetch_json(url).get("molecules") or []
    except FetchError:
        return None
    if not molecules:
        _MOL_CACHE[key] = None
        return None

    # Prefer an exact preferred-name match; otherwise the first hit.
    chosen = next(
        (m for m in molecules if (m.get("pref_name") or "").lower() == key),
        molecules[0],
    )
    hierarchy = chosen.get("molecule_hierarchy") or {}
    parent_id = hierarchy.get("parent_chembl_id") or chosen.get("molecule_chembl_id")
    pref = chosen.get("pref_name")
    phase = _phase_num(chosen.get("max_phase"))

    # If we fell back to a salt/child, fetch the parent record for clean fields.
    if parent_id and parent_id != chosen.get("molecule_chembl_id"):
        try:
            parent = fetch_json(f"{_BASE}/molecule/{parent_id}?format=json")
            pref = parent.get("pref_name") or pref
            phase = _phase_num(parent.get("max_phase")) if phase is None else phase
        except FetchError:
            pass

    out = {"chembl_id": parent_id, "pref_name": pref, "max_phase": phase}
    _MOL_CACHE[key] = out
    return out


def _get_target(target_chembl_id: str) -> dict:
    if target_chembl_id in _TARGET_CACHE:
        return _TARGET_CACHE[target_chembl_id]
    try:
        t = fetch_json(f"{_BASE}/target/{target_chembl_id}?format=json")
    except FetchError:
        t = {}
    info = {
        "name": t.get("pref_name"),
        "accessions": [
            c.get("accession")
            for c in (t.get("target_components") or [])
            if c.get("accession")
        ],
        "organism": t.get("organism"),
    }
    _TARGET_CACHE[target_chembl_id] = info
    return info


def _best_activity(molecule_id: str, target_id: str):
    """Return the most potent measured Ki/IC50/Kd (nM) for a molecule-target pair."""
    try:
        url = (
            f"{_BASE}/activity?molecule_chembl_id={molecule_id}"
            f"&target_chembl_id={target_id}"
            "&standard_type__in=Ki,IC50,Kd,EC50&format=json&limit=100"
        )
        acts = fetch_json(url).get("activities") or []
    except FetchError:
        return None
    best = None
    for a in acts:
        if (a.get("standard_units") or "").lower() != "nm":
            continue
        try:
            val = float(a.get("standard_value"))
        except (TypeError, ValueError):
            continue
        if best is None or val < best["value_nM"]:
            best = {
                "type": a.get("standard_type"),
                "relation": a.get("standard_relation") or "=",
                "value_nM": val,
            }
    return best


def _tokens(text: str) -> set:
    return {
        w for w in re.split(r"[^a-z0-9]+", (text or "").lower())
        if w and w not in _STOP and len(w) > 2
    }


def _name_matches(target_name: str, protein_title: str) -> bool:
    a, b = _tokens(target_name), _tokens(protein_title)
    shared = a & b
    if not shared:
        return False
    # A single very distinctive family token (protease, trypsin, anhydrase,
    # reductase, polymerase, ...) is enough; otherwise need 2+ shared tokens.
    if any(len(t) >= 7 for t in shared):
        return True
    strong = {t for t in shared if len(t) >= 5}
    return len(shared) >= 2 and len(strong) >= 1


def pharmacology(name: str, uniprots=None, protein_title: str = "") -> dict | None:
    """Curated ChEMBL pharmacology for a chemical + match vs the loaded protein.

    Returns None when the chemical isn't in ChEMBL. Otherwise a dict with the
    molecule, mechanisms of action (with best measured potency per target), and
    a `match` describing how the drug's targets relate to the loaded protein.
    """
    mol = _resolve_molecule(name)
    if not mol or not mol.get("chembl_id"):
        return None
    mol_id = mol["chembl_id"]
    uniprots = set(uniprots or [])

    mechanisms = []
    try:
        url = f"{_BASE}/mechanism?molecule_chembl_id={mol_id}&format=json&limit=8"
        raw = fetch_json(url).get("mechanisms") or []
    except FetchError:
        raw = []

    match = {"level": "none", "target_name": None, "best_activity": None}

    seen_targets = set()
    for m in raw:
        tid = m.get("target_chembl_id")
        if not tid or tid in seen_targets:
            continue
        seen_targets.add(tid)
        tinfo = _get_target(tid)
        best = _best_activity(mol_id, tid)
        mech = {
            "moa": m.get("mechanism_of_action"),
            "target_chembl_id": tid,
            "target_name": tinfo.get("name"),
            "organism": tinfo.get("organism"),
            "best_activity": best,
        }
        mechanisms.append(mech)

        # Cross-check this mechanism target against the loaded protein.
        if match["level"] != "uniprot":
            if uniprots & set(tinfo.get("accessions") or []):
                match = {"level": "uniprot", "target_name": tinfo.get("name"),
                         "best_activity": best}
            elif match["level"] == "none" and _name_matches(tinfo.get("name"), protein_title):
                match = {"level": "name", "target_name": tinfo.get("name"),
                         "best_activity": best}

        if len(mechanisms) >= 4:
            break

    return {
        "chembl_id": mol_id,
        "pref_name": mol.get("pref_name"),
        "max_phase": mol.get("max_phase"),
        "development_status": _PHASE.get(mol.get("max_phase"), "unknown"),
        "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{mol_id}/",
        "mechanisms": mechanisms,
        "match": match,
    }
