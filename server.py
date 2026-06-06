"""AtomScope local web server (Python stdlib only).

Serves the single-page frontend and a small JSON API that drives
structure loading, atomic interaction profiling, and chemical lookup.

Run:
    python server.py            # http://127.0.0.1:8010
    python server.py --port 8000
"""

from __future__ import annotations

import argparse
import json
import os
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import datetime

from atomscope import __version__ as ATOMSCOPE_VERSION
from atomscope import (
    chembl,
    docking,
    interactions,
    pdbparse,
    pockets,
    pubchem,
    rcsb,
    report,
)
from atomscope.http_util import FetchError

WEB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")

# Bounded in-memory cache of parsed structures: pdb_id -> (text, Structure, meta)
_CACHE: dict[str, tuple] = {}
_CACHE_LOCK = threading.Lock()
_CACHE_MAX = 16

# Pocket detection is expensive; cache results per PDB id.
_POCKET_CACHE: dict[str, list] = {}
_POCKET_LOCK = threading.Lock()

_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".ico": "image/x-icon",
}


def _load_structure(pdb_id: str):
    pid = rcsb.normalize_pdb_id(pdb_id)
    with _CACHE_LOCK:
        cached = _CACHE.get(pid)
    if cached:
        return cached

    text = rcsb.fetch_structure(pid)
    structure = pdbparse.parse_pdb(text)
    try:
        meta = rcsb.fetch_entry_metadata(pid)
    except FetchError:
        meta = {"pdb_id": pid, "title": None}
    entry = (text, structure, meta)

    with _CACHE_LOCK:
        if len(_CACHE) >= _CACHE_MAX:
            _CACHE.pop(next(iter(_CACHE)))
        _CACHE[pid] = entry
    return entry


_UNIPROT_CACHE: dict[str, list] = {}


def _get_uniprots(pdb_id: str) -> list:
    pid = rcsb.normalize_pdb_id(pdb_id)
    if pid in _UNIPROT_CACHE:
        return _UNIPROT_CACHE[pid]
    try:
        accs = rcsb.fetch_uniprot_accessions(pid)
    except FetchError:
        accs = []
    _UNIPROT_CACHE[pid] = accs
    return accs


def _get_pockets(pdb_id: str) -> list:
    pid = rcsb.normalize_pdb_id(pdb_id)
    with _POCKET_LOCK:
        cached = _POCKET_CACHE.get(pid)
    if cached is not None:
        return cached
    _text, structure, _meta = _load_structure(pid)
    found = pockets.detect_pockets(structure)
    with _POCKET_LOCK:
        if len(_POCKET_CACHE) >= _CACHE_MAX:
            _POCKET_CACHE.pop(next(iter(_POCKET_CACHE)))
        _POCKET_CACHE[pid] = found
    return found


def _methods_block(meta, site_label, center, search, ligand=None):
    """Assemble a reproducibility/methods record for a docking or screening run.

    Follows the docking literature's reporting recommendations: tool version,
    receptor identity, box definition, grid spacing, scoring model, search
    settings + random seed, ligand source/flexibility, and prep assumptions.
    """
    block = {
        "tool": f"AtomScope v{ATOMSCOPE_VERSION}",
        "run_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"
        ),
        "receptor": {
            "pdb_id": meta.get("pdb_id"),
            "title": meta.get("title"),
            "method": meta.get("experimental_method"),
            "resolution_A": meta.get("resolution_A"),
        },
        "receptor_prep": (
            "Rigid receptor; protein heavy atoms only — waters, ions, cofactors "
            "and any existing ligands are excluded from the scoring grid. No "
            "explicit hydrogens; heavy-atom geometry at implicit pH ~7."
        ),
        "site": site_label,
        "box": {
            "center": [round(c, 2) for c in center],
            "edge_A": round(2 * docking.GRID_HALF, 1),
            "grid_spacing_A": docking.SPACING,
            "translation_search_A": docking.TRANS_HALF,
        },
        "scoring": (
            "AutoDock-style grid-map empirical score: steric (smoothed "
            "Lennard-Jones) + hydrogen-bond + hydrophobic channels, "
            "trilinear-interpolated. Relative units (lower = better) — NOT "
            "calibrated to kcal/mol."
        ),
        "search": {
            "algorithm": "Monte-Carlo rigid-body, simulated-annealing acceptance",
            "seeds": search.get("seeds"),
            "mc_steps": search.get("mc_steps"),
            "random_seed": search.get("random_seed"),
        },
        "interaction_cutoffs_A": {
            "hydrogen_bond": interactions.HB_MAX,
            "salt_bridge": interactions.SALT_MAX,
            "hydrophobic": interactions.HYDRO_MAX,
            "metal_coordination": interactions.METAL_MAX,
            "aromatic_centroid": interactions.ARO_CENTROID_MAX,
        },
        "disclaimer": (
            "Research-only. Predicted poses and scores are geometric/empirical "
            "heuristics from a single static structure; not affinities, not "
            "clinical guidance. Validate with orthogonal evidence."
        ),
    }
    if ligand is not None:
        block["ligand"] = {
            "source": "PubChem",
            "cid": ligand.get("cid"),
            "conformer": ligand.get("conformer"),
            "n_heavy_atoms": ligand.get("n_heavy_atoms"),
            "flexibility": "rigid (single PubChem 3D conformer)",
        }
    return block


def _resolve_dock_site(pdb_id, structure, comp_raw, pocket_raw):
    """Return (center, label, ref_component) for a component or detected pocket.

    ref_component is the crystallographic Component when docking into an existing
    ligand site (enables redock RMSD), else None. Raises ValueError on bad input.
    """
    if comp_raw != "":
        try:
            idx = int(comp_raw)
        except ValueError:
            raise ValueError("'comp' must be an integer index")
        if idx < 0 or idx >= len(structure.components):
            raise ValueError("Component index out of range")
        comp = structure.components[idx]
        return docking.component_center(comp), comp.label, comp
    if pocket_raw != "":
        try:
            pidx = int(pocket_raw)
        except ValueError:
            raise ValueError("'pocket' must be an integer index")
        found = _get_pockets(pdb_id)
        match = next((p for p in found if p["index"] == pidx), None)
        if match is None:
            raise ValueError("Pocket index out of range")
        return tuple(match["center"]), f"detected pocket #{pidx} ({match['volume_A3']} Å³)", None
    raise ValueError("Need a docking site ('comp' or 'pocket')")


def _components_json(structure) -> list[dict]:
    out = []
    for i, c in enumerate(structure.components):
        out.append(
            {
                "index": i,
                "label": c.label,
                "res_name": c.res_name,
                "chain": c.chain,
                "res_seq": c.res_seq,
                "kind": c.kind,
                "atom_count": len(c.atoms),
            }
        )
    return out


class Handler(BaseHTTPRequestHandler):
    server_version = "AtomScope/0.1"

    def log_message(self, fmt, *args):  # quieter console
        pass

    # ---- helpers ------------------------------------------------------
    def _send_json(self, payload, status=200):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_error_json(self, message, status=400):
        self._send_json({"error": message}, status=status)

    def _serve_static(self, path):
        if path in ("/", ""):
            path = "/index.html"
        rel = path.lstrip("/")
        full = os.path.normpath(os.path.join(WEB_DIR, rel))
        if not full.startswith(WEB_DIR) or not os.path.isfile(full):
            self.send_error(404, "Not found")
            return
        ext = os.path.splitext(full)[1].lower()
        ctype = _CONTENT_TYPES.get(ext, "application/octet-stream")
        with open(full, "rb") as fh:
            body = fh.read()
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    # ---- routing ------------------------------------------------------
    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path
        qs = parse_qs(parsed.query)

        try:
            if path == "/api/analyze":
                return self._api_analyze(qs)
            if path == "/api/interactions":
                return self._api_interactions(qs)
            if path == "/api/chemical":
                return self._api_chemical(qs)
            if path == "/api/pockets":
                return self._api_pockets(qs)
            if path == "/api/dock":
                return self._api_dock(qs)
            if path == "/api/screen":
                return self._api_screen(qs)
            if path == "/api/search":
                return self._api_search(qs)
            return self._serve_static(path)
        except FetchError as exc:
            return self._send_error_json(str(exc), status=502)
        except Exception as exc:  # noqa: BLE001 - surface as JSON for the UI
            return self._send_error_json(f"Internal error: {exc}", status=500)

    # ---- API endpoints ------------------------------------------------
    def _api_analyze(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        text, structure, meta = _load_structure(pdb_id)
        return self._send_json(
            {
                "metadata": meta,
                "chains": structure.chains,
                "protein_atom_count": len(structure.protein_atoms),
                "components": _components_json(structure),
                "pdb_data": text,
            }
        )

    def _api_interactions(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        idx_raw = (qs.get("comp") or [""])[0]
        if not pdb_id or idx_raw == "":
            return self._send_error_json("Missing 'pdb' or 'comp' parameter")
        try:
            idx = int(idx_raw)
        except ValueError:
            return self._send_error_json("'comp' must be an integer index")

        _text, structure, meta = _load_structure(pdb_id)
        if idx < 0 or idx >= len(structure.components):
            return self._send_error_json("Component index out of range", status=404)
        component = structure.components[idx]
        profile = interactions.profile_component(structure, component)
        summary = report.summarize(profile, meta)
        return self._send_json({"profile": profile, "report": summary})

    def _api_chemical(self, qs):
        query = (qs.get("q") or [""])[0].strip()
        if not query:
            return self._send_error_json("Missing 'q' parameter")
        pdb_id = (qs.get("pdb") or [""])[0].strip()
        compound = pubchem.lookup_compound(query)

        # With a loaded protein, cross-reference ChEMBL pharmacology + measured
        # activity against that target; otherwise just the basic drug status.
        if pdb_id:
            try:
                _t, _s, meta = _load_structure(pdb_id)
                uniprots = _get_uniprots(pdb_id)
                compound["pharmacology"] = chembl.pharmacology(
                    query, uniprots, meta.get("title") or ""
                )
            except FetchError:
                compound["pharmacology"] = None
        else:
            compound["chembl"] = chembl.lookup_molecule(query)
        return self._send_json(compound)

    def _api_dock(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        chem = (qs.get("chem") or [""])[0].strip()
        comp_raw = (qs.get("comp") or [""])[0]
        pocket_raw = (qs.get("pocket") or [""])[0]
        if not pdb_id or not chem or (comp_raw == "" and pocket_raw == ""):
            return self._send_error_json(
                "Need 'pdb', 'chem', and a site ('comp' or 'pocket')"
            )

        _text, structure, meta = _load_structure(pdb_id)

        try:
            center, site_label, ref_component = _resolve_dock_site(
                pdb_id, structure, comp_raw, pocket_raw
            )
        except ValueError as exc:
            return self._send_error_json(str(exc))

        # Resolve the chemical and fetch a real 3D conformer from PubChem.
        compound = pubchem.lookup_compound(chem)
        if not compound.get("cid"):
            return self._send_error_json(f"Could not resolve chemical '{chem}'")
        lig = pubchem.fetch_3d_atoms(compound["cid"])

        res_name = (compound.get("molecular_formula") or "LIG")[:3].upper()

        pose = docking.dock(structure, lig["atoms"], center)
        comp = docking.pose_to_component(pose, res_name)
        profile = interactions.profile_component(structure, comp)
        summary = report.summarize(profile, meta)

        # Redock RMSD: only meaningful when the same molecule is crystallised here.
        redock_rmsd = None
        if ref_component is not None:
            ref_heavy = [a for a in ref_component.atoms if a.element != "H"]
            if len(ref_heavy) == len(lig["atoms"]):
                redock_rmsd = docking.rmsd_to_reference(pose, ref_component)

        # ChEMBL pharmacology + measured activity vs this target (validation).
        try:
            pharmacology = chembl.pharmacology(
                chem, _get_uniprots(pdb_id), meta.get("title") or ""
            )
        except FetchError:
            pharmacology = None

        return self._send_json(
            {
                "chemical": {
                    "cid": compound["cid"],
                    "name": compound.get("iupac_name") or chem,
                    "formula": compound.get("molecular_formula"),
                    "coord_source": lig["source"],
                    "n_heavy_atoms": len(lig["atoms"]),
                },
                "pocket": {"label": site_label, "center": pose["center"]},
                "docking": {
                    "score": pose["score"],
                    "ligand_efficiency": pose["ligand_efficiency"],
                    "box_half": pose["box_half"],
                    "search": pose["search"],
                    "redock_rmsd": redock_rmsd,
                },
                "pose_pdb": docking.pose_to_pdb(pose, res_name),
                "profile": profile,
                "report": summary,
                "pharmacology": pharmacology,
                "methods": _methods_block(
                    meta,
                    site_label,
                    pose["center"],
                    pose["search"],
                    ligand={
                        "cid": compound.get("cid"),
                        "conformer": lig["source"],
                        "n_heavy_atoms": len(lig["atoms"]),
                    },
                ),
            }
        )

    def _api_screen(self, qs):
        import re

        pdb_id = (qs.get("pdb") or [""])[0]
        chems_raw = (qs.get("chems") or [""])[0]
        comp_raw = (qs.get("comp") or [""])[0]
        pocket_raw = (qs.get("pocket") or [""])[0]
        if not pdb_id or not chems_raw:
            return self._send_error_json("Need 'pdb' and 'chems' (comma-separated)")

        tokens = [t.strip() for t in re.split(r"[,;\n]+", chems_raw) if t.strip()]
        # De-duplicate, preserve order, cap to keep runtime bounded.
        seen = set()
        chem_list = []
        for t in tokens:
            if t.lower() not in seen:
                seen.add(t.lower())
                chem_list.append(t)
        chem_list = chem_list[:10]
        if not chem_list:
            return self._send_error_json("No chemicals parsed from 'chems'")

        _text, structure, _meta = _load_structure(pdb_id)
        try:
            center, site_label, _ref = _resolve_dock_site(
                pdb_id, structure, comp_raw, pocket_raw
            )
        except ValueError as exc:
            return self._send_error_json(str(exc))

        # Build the scoring grid once, dock every ligand against it.
        grid = docking.build_grid(structure, center)
        results = []
        for token in chem_list:
            try:
                compound = pubchem.lookup_compound(token)
                if not compound.get("cid"):
                    raise FetchError(f"not found in PubChem")
                lig = pubchem.fetch_3d_atoms(compound["cid"])
                pose = docking.dock_with_grid(grid, lig["atoms"], center, seeds=160)
                res_name = (compound.get("molecular_formula") or "LIG")[:3].upper()
                comp = docking.pose_to_component(pose, res_name)
                profile = interactions.profile_component(structure, comp)
                top = [
                    f"{r['res_name']}{r['res_seq']}"
                    for r in profile["contact_residues"][:3]
                ]
                results.append(
                    {
                        "query": token,
                        "cid": compound["cid"],
                        "name": compound.get("iupac_name") or token,
                        "formula": compound.get("molecular_formula"),
                        "score": pose["score"],
                        "ligand_efficiency": pose["ligand_efficiency"],
                        "n_heavy_atoms": pose["n_heavy_atoms"],
                        "counts": profile["counts"],
                        "interaction_total": profile["interaction_total"],
                        "contact_residue_count": profile["contact_residue_count"],
                        "top_residues": top,
                        "coord_source": lig["source"],
                    }
                )
            except (FetchError, ValueError) as exc:
                results.append({"query": token, "error": str(exc)})

        # Rank: best (lowest) score first; failures last.
        results.sort(key=lambda r: r.get("score", float("inf")) if "error" not in r else float("inf"))
        for rank, r in enumerate(results, start=1):
            if "error" not in r:
                r["rank"] = rank

        methods = _methods_block(
            _meta,
            site_label,
            center,
            {"seeds": 160, "mc_steps": 40, "random_seed": 0},
        )
        return self._send_json(
            {
                "pdb_id": rcsb.normalize_pdb_id(pdb_id),
                "site": site_label,
                "count": len(results),
                "results": results,
                "methods": methods,
            }
        )

    def _api_pockets(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        found = _get_pockets(pdb_id)
        return self._send_json({"pockets": found, "count": len(found)})

    def _api_search(self, qs):
        query = (qs.get("q") or [""])[0].strip()
        if not query:
            return self._send_error_json("Missing 'q' parameter")
        return self._send_json({"results": rcsb.search_by_name(query, limit=10)})


def main():
    # Cloud hosts (Render/Railway/Fly/etc.) inject the port via $PORT and need
    # the server bound to all interfaces. Locally these default to 8010/127.0.0.1.
    env_port = os.environ.get("PORT")
    default_port = int(env_port) if env_port and env_port.isdigit() else 8010
    default_host = "0.0.0.0" if env_port else "127.0.0.1"

    parser = argparse.ArgumentParser(description="AtomScope web server")
    parser.add_argument("--port", type=int, default=default_port)
    parser.add_argument("--host", default=default_host)
    args = parser.parse_args()

    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"AtomScope running at http://{args.host}:{args.port}")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        httpd.shutdown()


if __name__ == "__main__":
    main()
