"""PubChem PUG-REST client: chemical lookup + druglikeness.

Handles drugs, chemicals, and single elements/ions by name or CID.
"""

from __future__ import annotations

import urllib.parse

from .http_util import FetchError, RateLimitError, fetch_json, fetch_text

_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Requested in priority order; SMILES field names have shifted across PubChem
# versions, so we degrade gracefully if a property set is rejected.
_PROP_SETS = [
    "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,"
    "HBondAcceptorCount,RotatableBondCount,Charge,IUPACName,SMILES",
    "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,"
    "HBondAcceptorCount,RotatableBondCount,Charge,IUPACName,ConnectivitySMILES",
    "MolecularFormula,MolecularWeight,XLogP,HBondDonorCount,"
    "HBondAcceptorCount,RotatableBondCount,IUPACName",
]


def _resolve_namespace(identifier: str) -> tuple[str, str]:
    ident = (identifier or "").strip()
    if not ident:
        raise FetchError("Empty chemical identifier")
    if ident.isdigit():
        return "cid", ident
    return "name", urllib.parse.quote(ident)


def lookup_compound(identifier: str) -> dict:
    """Look up a compound by name or CID and return properties + druglikeness."""
    namespace, value = _resolve_namespace(identifier)

    props = None
    last_err: Exception | None = None
    for prop_set in _PROP_SETS:
        url = f"{_BASE}/compound/{namespace}/{value}/property/{prop_set}/JSON"
        try:
            data = fetch_json(url)
            table = data.get("PropertyTable", {}).get("Properties", [])
            if table:
                props = table[0]
                break
        except RateLimitError:
            raise  # surface throttling clearly instead of "not found"
        except FetchError as exc:
            last_err = exc
            continue

    if props is None:
        raise FetchError(
            f"No PubChem compound found for '{identifier}'"
            + (f" ({last_err})" if last_err else "")
        )

    cid = props.get("CID")
    smiles = props.get("SMILES") or props.get("ConnectivitySMILES")

    result = {
        "cid": cid,
        "query": identifier,
        "iupac_name": props.get("IUPACName"),
        "molecular_formula": props.get("MolecularFormula"),
        "molecular_weight": _to_float(props.get("MolecularWeight")),
        "xlogp": _to_float(props.get("XLogP")),
        "tpsa": _to_float(props.get("TPSA")),
        "h_bond_donors": props.get("HBondDonorCount"),
        "h_bond_acceptors": props.get("HBondAcceptorCount"),
        "rotatable_bonds": props.get("RotatableBondCount"),
        "formal_charge": props.get("Charge"),
        "smiles": smiles,
        "image_url": f"{_BASE}/compound/cid/{cid}/PNG" if cid else None,
        "pubchem_url": (
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else None
        ),
    }
    result["druglikeness"] = _lipinski(result)
    result["rules"] = _extended_rules(result)
    return result


def fetch_3d_atoms(cid) -> dict:
    """Fetch real 3D coordinates for a CID from PubChem (heavy atoms only).

    Returns {"atoms": [{element, x, y, z}], "source": "3d"|"2d", "cid": cid}.
    Falls back to 2D coordinates (flagged) if no 3D conformer exists.
    """
    cid = str(cid)
    for record_type in ("3d", "2d"):
        url = f"{_BASE}/compound/cid/{cid}/SDF?record_type={record_type}"
        try:
            sdf = fetch_text(url)
        except RateLimitError:
            raise  # surface throttling clearly
        except FetchError:
            continue
        atoms = _parse_sdf_atoms(sdf)
        if atoms:
            return {"atoms": atoms, "source": record_type, "cid": int(cid)}
    raise FetchError(f"No 3D/2D coordinates available for CID {cid}")


def _parse_sdf_atoms(sdf: str) -> list[dict]:
    lines = sdf.splitlines()
    if len(lines) < 4:
        return []
    counts = lines[3]
    try:
        n_atoms = int(counts[0:3])
    except ValueError:
        return []
    atoms = []
    for i in range(4, 4 + n_atoms):
        if i >= len(lines):
            break
        line = lines[i]
        try:
            x = float(line[0:10])
            y = float(line[10:20])
            z = float(line[20:30])
            element = line[31:34].strip().upper()
        except ValueError:
            continue
        if element == "H":
            continue
        atoms.append({"element": element, "x": x, "y": y, "z": z})
    return atoms


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _extended_rules(p: dict) -> dict:
    """Druglikeness rules beyond Lipinski, computed from PubChem descriptors.

    Veber (oral bioavailability), Egan (absorption), lead-likeness, and a
    BOILED-Egg-style absorption/BBB read (XLogP used as a WLOGP proxy). All are
    transparent rule-of-thumb filters, not trained models.
    """
    mw = p.get("molecular_weight")
    logp = p.get("xlogp")
    tpsa = p.get("tpsa")
    rot = p.get("rotatable_bonds")

    def rule(passed, criteria):
        return {"pass": (None if passed is None else bool(passed)), "criteria": criteria}

    veber = None if (rot is None or tpsa is None) else (rot <= 10 and tpsa <= 140)
    egan = None if (tpsa is None or logp is None) else (tpsa <= 131.6 and logp <= 5.88)
    lead = (
        None
        if (mw is None or logp is None or rot is None)
        else (mw <= 350 and logp <= 3.5 and rot <= 7)
    )

    gi = bbb = None
    if tpsa is not None and logp is not None:
        gi = "High" if (tpsa <= 130 and -0.5 <= logp <= 6.0) else "Low"
        bbb = "Yes" if (tpsa <= 79 and 0.4 <= logp <= 6.0) else "No"

    return {
        "veber": rule(veber, "rotatable ≤10 and TPSA ≤140"),
        "egan": rule(egan, "TPSA ≤131.6 and XLogP ≤5.88"),
        "lead_like": rule(lead, "MW ≤350, XLogP ≤3.5, rotatable ≤7"),
        "absorption": {
            "gi_absorption": gi,
            "bbb_permeant": bbb,
            "model": "BOILED-Egg style (approx; XLogP as WLogP proxy)",
        },
    }


def _lipinski(props: dict) -> dict:
    """Lipinski Rule of Five evaluation (research heuristic only)."""
    violations = []
    mw = props.get("molecular_weight")
    logp = props.get("xlogp")
    hbd = props.get("h_bond_donors")
    hba = props.get("h_bond_acceptors")

    if mw is not None and mw > 500:
        violations.append(f"MW {mw:.0f} > 500")
    if logp is not None and logp > 5:
        violations.append(f"XLogP {logp:.1f} > 5")
    if hbd is not None and hbd > 5:
        violations.append(f"H-bond donors {hbd} > 5")
    if hba is not None and hba > 10:
        violations.append(f"H-bond acceptors {hba} > 10")

    return {
        "rule": "Lipinski Rule of Five",
        "violations": violations,
        "violation_count": len(violations),
        "drug_like": len(violations) <= 1,
    }
