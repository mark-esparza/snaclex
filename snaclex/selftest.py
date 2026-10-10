"""Upstream connectivity self-test: does every integrated database actually work?

SnaCleX has no database of its own. Every result it produces is assembled live
from five public services, and a failure in any one of them degrades or breaks a
feature without necessarily producing an obvious error:

    RCSB PDB      files.rcsb.org   structure coordinates (PDB + mmCIF)
                  data.rcsb.org    entry metadata, UniProt mapping (GraphQL)
                  search.rcsb.org  full-text entry search
    PubChem       pubchem.ncbi.nlm.nih.gov   compound properties, 3D conformers
    ChEMBL        www.ebi.ac.uk/chembl       pharmacology and measured activity
    InterPro/Pfam www.ebi.ac.uk/interpro     family alignments for conservation

This module checks all of them, and checks them *properly*. A reachability ping
only proves a host answered; it does not prove the integration still works. Each
check therefore calls SnaCleX's own client for that source and asserts a known
fact about the result — aspirin is CID 2244, crambin has a few hundred protein
atoms — so a silently changed response shape fails here rather than in a user's
analysis.

Run it against any deployment::

    python -m snaclex.selftest            # human-readable table
    python -m snaclex.selftest --json     # machine-readable
    curl https://<host>/api/selftest      # the same checks, on the live host

Exit status is 0 only when every check passes, so it can gate a deploy.
"""

from __future__ import annotations

import datetime
import json
import time

# Fixed, tiny, long-lived reference entries. Crambin (46 residues) is one of the
# smallest protein structures in the PDB, and aspirin is about as stable an
# identifier as PubChem/ChEMBL have, so these checks stay cheap and durable.
REF_PDB = "1CRN"
REF_UNIPROT = "P01542"      # crambin
REF_CHEM = "aspirin"
REF_CID = 2244


def _check_rcsb_structure():
    from . import pdbparse, rcsb

    text = rcsb.fetch_structure(REF_PDB)
    structure = pdbparse.parse_structure(text)
    n = len(structure.protein_atoms)
    if n < 100:
        raise ValueError(f"{REF_PDB} parsed to only {n} protein atoms")
    return f"{REF_PDB}: {n} protein atoms parsed"


def _check_rcsb_metadata():
    from . import rcsb

    meta = rcsb.fetch_entry_metadata(REF_PDB)
    if not meta.get("title"):
        raise ValueError("entry metadata has no title")
    return f"{REF_PDB}: {meta['title'][:60]}"


def _check_rcsb_search():
    from . import rcsb

    hits = rcsb.search_by_name("crambin", limit=3)
    if not hits:
        raise ValueError("full-text search returned no hits for 'crambin'")
    return f"{len(hits)} hit(s); first = {hits[0].get('pdb_id') or hits[0]}"


def _check_rcsb_uniprot():
    from . import rcsb

    accs = rcsb.fetch_uniprot_accessions(REF_PDB)
    if not accs:
        raise ValueError(f"no UniProt accession mapped for {REF_PDB}")
    return f"{REF_PDB} -> {', '.join(accs)}"


def _check_pubchem_properties():
    from . import pubchem

    c = pubchem.lookup_compound(REF_CHEM)
    if c.get("cid") != REF_CID:
        raise ValueError(f"expected CID {REF_CID} for {REF_CHEM}, got {c.get('cid')}")
    if not c.get("molecular_formula"):
        raise ValueError("compound record has no molecular formula")
    return f"{REF_CHEM} = CID {c['cid']} ({c['molecular_formula']})"


def _check_pubchem_conformer():
    from . import pubchem

    lig = pubchem.fetch_3d_atoms(REF_CID)
    n = len(lig.get("atoms") or [])
    if n < 5:
        raise ValueError(f"conformer for CID {REF_CID} has only {n} atoms")
    return f"CID {REF_CID}: {n} heavy atoms ({lig.get('source')} conformer)"


def _check_chembl():
    from . import chembl

    mol = chembl.lookup_molecule(REF_CHEM)
    if not mol or not mol.get("chembl_id"):
        raise ValueError(f"no ChEMBL molecule resolved for {REF_CHEM}")
    return f"{REF_CHEM} = {mol['chembl_id']}"


def _check_interpro():
    from . import evolution

    fam = evolution.pfam_for_uniprot(REF_UNIPROT)
    if not fam or not fam.get("pfam"):
        # pfam_for_uniprot swallows FetchError and returns None, so an
        # unreachable host and a genuinely family-less accession look alike
        # from here. Say so rather than asserting which one it was.
        raise ValueError(
            f"no Pfam family for {REF_UNIPROT} — host unreachable, or the "
            f"InterPro response shape changed"
        )
    return f"{REF_UNIPROT} -> {fam['pfam']} ({fam.get('name') or '?'})"


# (id, source, host, what it proves, callable)
CHECKS = [
    ("rcsb_structure", "RCSB PDB", "files.rcsb.org",
     "structure coordinates download + parse", _check_rcsb_structure),
    ("rcsb_metadata", "RCSB PDB", "data.rcsb.org",
     "entry metadata (title, method, resolution)", _check_rcsb_metadata),
    ("rcsb_search", "RCSB PDB", "search.rcsb.org",
     "full-text entry search", _check_rcsb_search),
    ("rcsb_uniprot", "RCSB PDB", "data.rcsb.org",
     "UniProt mapping via GraphQL (feeds conservation)", _check_rcsb_uniprot),
    ("pubchem_properties", "PubChem", "pubchem.ncbi.nlm.nih.gov",
     "compound properties + druglikeness", _check_pubchem_properties),
    ("pubchem_conformer", "PubChem", "pubchem.ncbi.nlm.nih.gov",
     "3D conformer download (feeds docking)", _check_pubchem_conformer),
    ("chembl", "ChEMBL", "www.ebi.ac.uk",
     "molecule lookup (feeds pharmacology + measured activity)", _check_chembl),
    ("interpro", "InterPro / Pfam", "www.ebi.ac.uk",
     "family lookup (feeds conservation scoring)", _check_interpro),
]


def run_checks(checks=None, timer=time.monotonic) -> dict:
    """Run every upstream check and return a structured report.

    Never raises: a failing upstream is a *result*, not an error, since the
    whole point is to report which integrations are down.
    """
    from .http_util import fail_fast

    results = []
    for check_id, source, host, proves, fn in (checks or CHECKS):
        started = timer()
        entry = {
            "id": check_id,
            "source": source,
            "host": host,
            "proves": proves,
        }
        try:
            # A diagnostic reports a down host quickly; it does not sit through
            # the retry budget a user-facing analysis is entitled to.
            with fail_fast():
                entry["detail"] = fn()
            entry["ok"] = True
        except Exception as exc:  # noqa: BLE001 - a down upstream is the finding
            entry["ok"] = False
            entry["error"] = f"{type(exc).__name__}: {exc}"
        entry["ms"] = int((timer() - started) * 1000)
        results.append(entry)

    passed = [r for r in results if r["ok"]]
    failed = [r for r in results if not r["ok"]]
    by_source: dict[str, dict] = {}
    for r in results:
        agg = by_source.setdefault(r["source"], {"ok": 0, "failed": 0})
        agg["ok" if r["ok"] else "failed"] += 1

    return {
        "tool": "SnaCleX upstream self-test",
        "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"
        ),
        "n_checks": len(results),
        "n_passed": len(passed),
        "n_failed": len(failed),
        "all_ok": not failed,
        "by_source": by_source,
        "checks": results,
        "note": (
            "Each check calls SnaCleX's own client for that source and asserts a "
            "known fact about the response, so it fails on a changed response "
            "shape as well as on an unreachable host."
        ),
    }


def format_report(report) -> str:
    lines = [
        f"{report['tool']} — {report['run_utc']}",
        f"{report['n_passed']}/{report['n_checks']} checks passed",
        "",
    ]
    for r in report["checks"]:
        mark = "ok  " if r["ok"] else "FAIL"
        lines.append(f"  [{mark}] {r['source']:<16} {r['host']:<26} {r['ms']:>6} ms")
        lines.append(f"         {r['proves']}")
        lines.append(f"         {r.get('detail') or r.get('error')}")
    if not report["all_ok"]:
        lines += [
            "",
            "A failing check means that feature is degraded or broken on this host:",
            "  RCSB          -> no structures can be loaded at all",
            "  PubChem       -> no chemical lookup, no docking (needs a conformer)",
            "  ChEMBL        -> no pharmacology or measured-activity cross-reference",
            "  InterPro/Pfam -> no conservation scoring",
            "",
            "If the host has no outbound internet, these will all fail. Check the",
            "platform's egress rules before looking for a bug in SnaCleX.",
        ]
    return "\n".join(lines)


def _main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(prog="python -m snaclex.selftest")
    parser.add_argument("--json", action="store_true", help="emit JSON")
    args = parser.parse_args(argv)

    report = run_checks()
    print(json.dumps(report, indent=2) if args.json else format_report(report))
    return 0 if report["all_ok"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
