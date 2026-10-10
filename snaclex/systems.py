"""Curated biological systems: named panels of targets and ligands.

A *system* is a reusable, versioned description of one piece of biology -- the
enzymes, receptors and transporters that act on a related set of small
molecules -- so a panel can be run from a single id instead of being retyped
structure by structure. ``data/systems/*.json`` holds the curations;
:func:`load_system` validates and loads them.

Verification is a first-class field, not a footnote
---------------------------------------------------
A curated system asserts facts that a reader cannot check by running the code:
that PDB 1ABC really is the enzyme named, that it really contains the ligand
whose site is being docked into, that the UniProt accession matches. Those
assertions are exactly the kind that got flagged as unverifiable.

So every target carries ``verified: false`` until a real RCSB round-trip stamps
it, and :func:`load_system` surfaces an aggregate ``verification`` block that
callers are expected to display. Run::

    python -m snaclex.systems verify catecholamine --write

on a machine with network access. That fetches each entry, confirms the title,
records resolution/method, and checks the named site ligand is actually present
in the deposited coordinates -- then writes the evidence back into the JSON.
An unverified system is a *draft*, and the API says so on every response.
"""

from __future__ import annotations

import datetime
import json
import os

_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "systems")

# A site is addressed by the resname of a co-crystallized ligand ("dock where
# this known ligand sits") or by rank in SnaCleX's own pocket detection
# ("dock into the Nth detected cavity"). Resname addressing is preferred: it is
# stable across SnaCleX versions, whereas pocket ranks are not.
_SITE_KEYS = ("ligand", "pocket")


class SystemError_(ValueError):
    """Malformed or unknown system definition."""


def _require(doc, key, where):
    if key not in doc or doc[key] in (None, ""):
        raise SystemError_(f"{where}: missing required field '{key}'")
    return doc[key]


def validate_system(doc) -> dict:
    """Validate a system document's shape; return it unchanged.

    Raises :class:`SystemError_` with a specific message on the first problem,
    so a bad curation fails loudly at load rather than mid-panel.
    """
    if not isinstance(doc, dict):
        raise SystemError_("system must be a JSON object")
    sid = _require(doc, "id", "system")
    _require(doc, "name", f"system '{sid}'")

    targets = doc.get("targets")
    if not isinstance(targets, list) or not targets:
        raise SystemError_(f"system '{sid}': 'targets' must be a non-empty list")
    ligands = doc.get("ligands")
    if not isinstance(ligands, list) or not ligands:
        raise SystemError_(f"system '{sid}': 'ligands' must be a non-empty list")

    seen = set()
    for target in targets:
        where = f"system '{sid}' target"
        tid = _require(target, "id", where)
        if tid in seen:
            raise SystemError_(f"system '{sid}': duplicate target id '{tid}'")
        seen.add(tid)
        _require(target, "pdb", f"{where} '{tid}'")
        site = target.get("site")
        if not isinstance(site, dict) or not any(k in site for k in _SITE_KEYS):
            raise SystemError_(
                f"{where} '{tid}': 'site' must be an object with one of "
                f"{_SITE_KEYS}"
            )

    seen = set()
    for ligand in ligands:
        where = f"system '{sid}' ligand"
        lid = _require(ligand, "id", where)
        if lid in seen:
            raise SystemError_(f"system '{sid}': duplicate ligand id '{lid}'")
        seen.add(lid)
        _require(ligand, "query", f"{where} '{lid}'")

    return doc


def verification_summary(doc) -> dict:
    """Aggregate per-target verification state into a displayable block."""
    targets = doc.get("targets") or []
    verified = [t for t in targets if t.get("verified") is True]
    return {
        "n_targets": len(targets),
        "n_verified": len(verified),
        "fully_verified": len(verified) == len(targets) and bool(targets),
        "last_verified_utc": doc.get("last_verified_utc"),
        "unverified_target_ids": [
            t.get("id") for t in targets if t.get("verified") is not True
        ],
        "note": (
            "Every target's PDB entry has been checked against RCSB."
            if len(verified) == len(targets) and targets else
            "DRAFT CURATION: one or more structures have not been checked "
            "against RCSB. Structure identity, the named site ligand and the "
            "UniProt mapping are unconfirmed assertions by the curator. Run "
            "`python -m snaclex.systems verify <id> --write` before relying "
            "on, citing or publishing results from this system."
        ),
    }


def _path(system_id) -> str:
    safe = str(system_id).strip().lower()
    # Defend the data dir against traversal: ids are flat slugs.
    if not safe or not all(ch.isalnum() or ch in "-_" for ch in safe):
        raise SystemError_(f"Invalid system id '{system_id}'")
    return os.path.join(_DATA_DIR, f"{safe}.json")


def load_system(system_id) -> dict:
    """Load, validate and annotate one curated system."""
    path = _path(system_id)
    try:
        with open(path, encoding="utf-8") as fh:
            doc = json.load(fh)
    except FileNotFoundError:
        raise SystemError_(f"Unknown system '{system_id}'") from None
    except ValueError as exc:
        raise SystemError_(f"System '{system_id}' is not valid JSON: {exc}") from None
    validate_system(doc)
    doc["verification"] = verification_summary(doc)
    return doc


def list_systems() -> list[dict]:
    """Return a compact catalog of every curated system that loads cleanly."""
    out = []
    try:
        names = sorted(os.listdir(_DATA_DIR))
    except OSError:
        return out
    for name in names:
        if not name.endswith(".json"):
            continue
        try:
            doc = load_system(name[: -len(".json")])
        except SystemError_:
            continue  # a broken curation shouldn't break the catalog
        out.append({
            "id": doc["id"],
            "name": doc["name"],
            "description": doc.get("description"),
            "version": doc.get("version"),
            "n_targets": len(doc["targets"]),
            "n_ligands": len(doc["ligands"]),
            "verification": doc["verification"],
        })
    return out


def panel_specs(doc, target_ids=None, ligand_ids=None):
    """Project a system document into (targets, ligands) specs for `panel`.

    Optional id filters let a caller run a slice of a large system.
    """
    targets = doc["targets"]
    ligands = doc["ligands"]
    if target_ids:
        wanted = set(target_ids)
        targets = [t for t in targets if t["id"] in wanted]
        if not targets:
            raise SystemError_("No targets matched the requested target_ids")
    if ligand_ids:
        wanted = set(ligand_ids)
        ligands = [lig for lig in ligands if lig["id"] in wanted]
        if not ligands:
            raise SystemError_("No ligands matched the requested ligand_ids")
    return targets, ligands


# ---------------------------------------------------------------------------
# Verification CLI (needs network; see module docstring)
# ---------------------------------------------------------------------------

def verify_system(system_id, *, fetch_metadata=None, fetch_structure=None,
                  parse_structure=None, fetch_uniprots=None, write=False) -> dict:
    """Check every target in a system against RCSB and report findings.

    Three independent assertions are checked per target: that the entry exists
    and what it is actually called, that the declared UniProt accession really
    maps to it, and (for ligand-addressed sites) that the named ligand is
    genuinely present in the deposited coordinates.

    The RCSB/parse functions are injected so this is testable offline. With
    ``write=True`` the per-target ``verified`` flags and evidence are persisted
    back into the system's JSON file.
    """
    from . import pdbparse, rcsb

    fetch_metadata = fetch_metadata or rcsb.fetch_entry_metadata
    fetch_structure = fetch_structure or rcsb.fetch_structure
    parse_structure = parse_structure or pdbparse.parse_structure
    fetch_uniprots = fetch_uniprots or rcsb.fetch_uniprot_accessions

    path = _path(system_id)
    with open(path, encoding="utf-8") as fh:
        doc = json.load(fh)
    validate_system(doc)

    findings = []
    for target in doc["targets"]:
        pdb = target["pdb"]
        finding = {"target_id": target["id"], "pdb": pdb, "problems": []}
        try:
            meta = fetch_metadata(pdb)
        except Exception as exc:  # noqa: BLE001 - report, don't abort the sweep
            finding["problems"].append(f"could not fetch entry: {exc}")
            target["verified"] = False
            findings.append(finding)
            continue

        finding["title"] = meta.get("title")
        finding["resolution_A"] = meta.get("resolution_A")
        finding["method"] = meta.get("experimental_method")

        # Confirm the declared UniProt accession really maps to this entry.
        declared = target.get("uniprot")
        accessions = None
        if declared:
            try:
                accessions = fetch_uniprots(pdb)
                finding["uniprot_accessions"] = accessions
                if accessions and declared not in accessions:
                    finding["problems"].append(
                        f"declared UniProt '{declared}' is not among {pdb}'s "
                        f"accessions ({', '.join(accessions)}) -- the entry may "
                        f"be a different protein or a non-human orthologue"
                    )
            except Exception as exc:  # noqa: BLE001
                finding["problems"].append(f"could not fetch UniProt mapping: {exc}")

        # Download coordinates once if either the site ligand or any declared
        # cofactor needs confirming against what was actually deposited.
        wanted = (target.get("site") or {}).get("ligand")
        declared_cofactors = [str(c).upper() for c in (target.get("cofactors") or [])]
        present = None
        cofactors_found = None
        if wanted or declared_cofactors:
            try:
                structure = parse_structure(fetch_structure(pdb))
                all_het = sorted({c.res_name.upper() for c in structure.components})
                present = sorted({
                    c.res_name.upper() for c in structure.ligand_components
                })
                finding["ligands_present"] = present

                if wanted and wanted.upper() not in present:
                    finding["problems"].append(
                        f"site ligand '{wanted}' not found in {pdb} "
                        f"(present: {', '.join(present) or 'none'})"
                    )

                if declared_cofactors:
                    cofactors_found = [c for c in declared_cofactors if c in all_het]
                    finding["cofactors_found"] = cofactors_found
                    # Alternatives are declared per target (SAM or SAH, PLP or
                    # LLP, ...), so none matching is the real problem.
                    if not cofactors_found:
                        finding["problems"].append(
                            f"none of the declared cofactors "
                            f"({', '.join(declared_cofactors)}) are present in "
                            f"{pdb}; this site would be scored as an empty "
                            f"cavity (deposited heterocomponents: "
                            f"{', '.join(all_het) or 'none'})"
                        )
            except Exception as exc:  # noqa: BLE001
                finding["problems"].append(f"could not parse coordinates: {exc}")

        ok = not finding["problems"]
        target["verified"] = ok
        target["verification"] = {
            "checked_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
                "%Y-%m-%d %H:%M UTC"
            ),
            "title": meta.get("title"),
            "resolution_A": meta.get("resolution_A"),
            "experimental_method": meta.get("experimental_method"),
            "uniprot_accessions": accessions,
            "ligands_present": present,
            "cofactors_found": cofactors_found,
            "problems": finding["problems"],
        }
        findings.append(finding)

    doc["last_verified_utc"] = datetime.datetime.now(
        datetime.timezone.utc
    ).strftime("%Y-%m-%d %H:%M UTC")

    if write:
        doc.pop("verification", None)  # derived; never persisted
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        os.replace(tmp, path)

    return {
        "system": doc["id"],
        "findings": findings,
        "n_ok": sum(1 for f in findings if not f["problems"]),
        "n_problem": sum(1 for f in findings if f["problems"]),
        "written": bool(write),
    }


def _main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(prog="python -m snaclex.systems")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list", help="list curated systems")
    show = sub.add_parser("show", help="print one system")
    show.add_argument("system")
    ver = sub.add_parser("verify", help="check a system's structures against RCSB")
    ver.add_argument("system")
    ver.add_argument("--write", action="store_true",
                     help="persist verification results into the JSON")
    args = parser.parse_args(argv)

    if args.cmd == "list":
        for entry in list_systems():
            flag = "ok " if entry["verification"]["fully_verified"] else "DRAFT"
            print(f"[{flag}] {entry['id']:<16} {entry['n_targets']:>2}T x "
                  f"{entry['n_ligands']:>2}L  {entry['name']}")
        return 0

    if args.cmd == "show":
        print(json.dumps(load_system(args.system), indent=2))
        return 0

    result = verify_system(args.system, write=args.write)
    for finding in result["findings"]:
        if finding["problems"]:
            print(f"FAIL {finding['target_id']} ({finding['pdb']})")
            for problem in finding["problems"]:
                print(f"       - {problem}")
        else:
            print(f"ok   {finding['target_id']} ({finding['pdb']}) "
                  f"{(finding.get('title') or '')[:70]}")
    print(f"\n{result['n_ok']} ok, {result['n_problem']} with problems"
          + (" (written)" if result["written"] else " (dry run; use --write)"))
    return 1 if result["n_problem"] else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
