"""SnaCleX local web server (Python stdlib only).

Copyright (c) 2026 Mark Esparza. All rights reserved. Proprietary — see LICENSE.


Serves the single-page frontend and a small JSON API that drives
structure loading, atomic interaction profiling, and chemical lookup.

Run:
    python server.py            # http://127.0.0.1:8010
    python server.py --port 8000
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import secrets
import hashlib
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

import datetime

from snaclex import __version__ as SNACLEX_VERSION
from snaclex import (
    antibody,
    apidocs,
    benchmark,
    bravo,
    cellxgene,
    chembl,
    docking,
    evolution,
    hla,
    interactions,
    jobs,
    models_esm,
    pdbparse,
    pockets,
    provenance,
    pubchem,
    rcsb,
    report,
    variants,
    vep,
)
from snaclex.http_util import FetchError

WEB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")

# ---------------------------------------------------------------------------
# Logging / observability
# ---------------------------------------------------------------------------
log = logging.getLogger("snaclex")

# Client IPs are logged only as a salted hash (the privacy policy says we don't
# track individuals). The salt is per-process unless pinned via SNACLEX_IP_SALT.
_IP_SALT = (os.environ.get("SNACLEX_IP_SALT") or secrets.token_hex(8)).encode()


def _ip_hash(ip: str) -> str:
    return hashlib.sha256(_IP_SALT + ip.encode("utf-8", "ignore")).hexdigest()[:12]


def configure_logging():
    level = os.environ.get("SNACLEX_LOG_LEVEL", "INFO").upper()
    logging.basicConfig(
        level=getattr(logging, level, logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )

# Bounded in-memory cache of parsed structures: pdb_id -> (text, Structure, meta)
_CACHE: dict[str, tuple] = {}
_CACHE_LOCK = threading.Lock()
_CACHE_MAX = 16

# Pocket detection is expensive; cache results per PDB id.
_POCKET_CACHE: dict[str, list] = {}
_POCKET_LOCK = threading.Lock()

# Docking scoring grids are expensive to build; cache per (pdb_id, site center)
# so repeat docks and every ligand in a screen reuse one grid.
_GRID_CACHE: dict[tuple, "docking.Grid"] = {}
_GRID_LOCK = threading.Lock()
_GRID_MAX = 16

# User-uploaded structures: id -> (text, Structure, meta). Kept only in memory,
# bounded, and never written to disk (see privacy policy).
_UPLOAD_CACHE: dict[str, tuple] = {}
_UPLOAD_LOCK = threading.Lock()
_UPLOAD_MAX = 8
MAX_UPLOAD_BYTES = 5_000_000


def _norm_id(pdb_id: str) -> str:
    """Normalize a structure id: pass through upload ids, else validate as PDB."""
    if pdb_id in _UPLOAD_CACHE:
        return pdb_id
    return rcsb.normalize_pdb_id(pdb_id)

_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".svg": "image/svg+xml",
    ".png": "image/png",
    ".ico": "image/x-icon",
}

# ---------------------------------------------------------------------------
# Security / abuse controls
# ---------------------------------------------------------------------------

# Content Security Policy tuned to exactly what the app loads: the 3Dmol.js
# viewer from its CDN (and the blob: web worker it spawns for surface meshes),
# PubChem 2D structure images, and same-origin everything else. Inline styles
# are allowed because the markup uses a few `style="…"` attributes.
_CSP = (
    "default-src 'self'; "
    "script-src 'self' https://3Dmol.org blob:; "
    "style-src 'self' 'unsafe-inline'; "
    "img-src 'self' data: blob: https://pubchem.ncbi.nlm.nih.gov; "
    "connect-src 'self'; "
    "worker-src 'self' blob:; "
    "font-src 'self'; object-src 'none'; base-uri 'self'; "
    "frame-ancestors 'none'; form-action 'self'"
)

# Compute-heavy GET endpoints get a tighter per-IP budget plus a global
# concurrency cap. Docking/screening are no longer here — they run through the
# async job queue (POST /api/jobs), which bounds concurrency via its worker pool.
EXPENSIVE_ENDPOINTS = {"/api/pockets", "/api/evolution", "/api/interface",
                       "/api/nucleic", "/api/hla", "/api/expression"}

MAX_QUERY_LEN = 200       # single chemical / search term
MAX_CHEMS_LEN = 2000      # batch-screening textarea

# A structure must fit PDB-format limits and stay interactive in the browser.
# Bigger entries (often mmCIF-only mega-assemblies) are rejected with a clear
# message instead of freezing the tab. Capped at the PDB serial limit (99,999).
try:
    MAX_STRUCTURE_ATOMS = min(99_999, int(os.environ.get("SNACLEX_MAX_ATOMS") or 99_999))
except ValueError:
    MAX_STRUCTURE_ATOMS = 99_999
MAX_STRUCTURE_BYTES = 30_000_000


def _env_float(name, default):
    raw = os.environ.get(name)
    try:
        return float(raw) if raw else default
    except ValueError:
        return default


class RateLimiter:
    """Thread-safe token bucket keyed by client identifier (IP).

    `capacity` tokens accrue at `rate` tokens/second; each allowed request
    spends `cost` tokens. Returns (allowed: bool, retry_after_seconds: float).
    """

    def __init__(self, rate, capacity, time_fn=time.monotonic, max_keys=4096):
        self.rate = rate
        self.capacity = capacity
        self._time = time_fn
        self._max_keys = max_keys
        self._lock = threading.Lock()
        self._buckets: dict[str, tuple] = {}

    def allow(self, key, cost=1.0):
        now = self._time()
        with self._lock:
            tokens, last = self._buckets.get(key, (self.capacity, now))
            tokens = min(self.capacity, tokens + (now - last) * self.rate)
            if tokens >= cost:
                self._buckets[key] = (tokens - cost, now)
                self._evict_if_needed()
                return True, 0.0
            self._buckets[key] = (tokens, now)
            retry = (cost - tokens) / self.rate if self.rate else 60.0
            return False, retry

    def _evict_if_needed(self):
        # Bound memory: drop the oldest-touched keys when the table grows large.
        if len(self._buckets) > self._max_keys:
            oldest = sorted(self._buckets, key=lambda k: self._buckets[k][1])
            for k in oldest[: self._max_keys // 4]:
                self._buckets.pop(k, None)


# General per-IP budget across all /api/* calls, plus a stricter one for the
# expensive compute endpoints. All tunable via env for ops.
_IP_LIMITER = RateLimiter(
    rate=_env_float("SNACLEX_RATE", 2.0),
    capacity=_env_float("SNACLEX_BURST", 60.0),
)
_HEAVY_LIMITER = RateLimiter(
    rate=_env_float("SNACLEX_HEAVY_RATE", 0.2),
    capacity=_env_float("SNACLEX_HEAVY_BURST", 6.0),
)
_HEAVY_SEM = threading.BoundedSemaphore(int(_env_float("SNACLEX_MAX_CONCURRENCY", 4)))

# Async job queue for long-running docking/screening (see snaclex/jobs.py).
JOBS = jobs.JobManager(
    max_workers=int(_env_float("SNACLEX_MAX_CONCURRENCY", 4)),
    ttl_seconds=int(_env_float("SNACLEX_JOB_TTL", 900)),
)


def clean_text(value, max_len=MAX_QUERY_LEN):
    """Strip NULs/control chars (keeping tab/newline) and cap a free-text param."""
    v = (value or "").replace("\x00", "")
    v = "".join(ch for ch in v if ch >= " " or ch in "\t\n")
    return v.strip()[:max_len]


def _check_structure_size(structure, label):
    """Reject structures too large for PDB format / interactive rendering."""
    n = len(structure.atoms)
    if n > MAX_STRUCTURE_ATOMS:
        raise FetchError(
            f"{label} has {n:,} atoms — too large for interactive analysis "
            f"(limit {MAX_STRUCTURE_ATOMS:,}). Try a single chain or the "
            f"asymmetric unit."
        )


def _load_full(pdb_id: str):
    """Fetch + parse a structure with no atom-count limit (cached).

    Used by analyze, which may then offer per-chain loading for big assemblies.
    Only a byte pre-cap guards against pathologically large downloads.
    """
    pid = rcsb.normalize_pdb_id(pdb_id)
    with _CACHE_LOCK:
        cached = _CACHE.get(pid)
    if cached:
        return cached

    text = rcsb.fetch_structure(pid)
    if len(text) > MAX_STRUCTURE_BYTES:
        raise FetchError(f"{pid} structure file is too large to load.")
    structure = pdbparse.parse_structure(text)
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


def _load_structure(pdb_id: str):
    """Return a size-OK, viewer-ready (PDB-format text) structure for analysis.

    Uploads and per-chain subsets live in the upload cache; plain PDB ids are
    fetched, size-checked, and re-serialized to PDB if they came as mmCIF.
    """
    with _UPLOAD_LOCK:
        up = _UPLOAD_CACHE.get(pdb_id)
    if up is not None:
        return up

    raw, structure, meta = _load_full(pdb_id)
    _check_structure_size(structure, meta.get("pdb_id") or pdb_id)
    text = pdbparse.to_pdb(structure) if "_atom_site." in raw else raw
    return (text, structure, meta)


def _cache_subset(parent_pid, chain, structure, meta):
    """Store a per-chain subset in the upload cache; return its synthetic id."""
    sub = pdbparse.subset_chain(structure, chain)
    _check_structure_size(sub, f"{parent_pid} chain {chain}")
    sub_id = f"{parent_pid}-{chain}"
    title = (meta.get("title") or parent_pid)
    sub_meta = {
        **meta,
        "pdb_id": sub_id,
        "source_pdb_id": parent_pid,
        "chain": chain,
        "title": f"{title} — chain {chain}",
        "subset": True,
    }
    entry = (pdbparse.to_pdb(sub), sub, sub_meta)
    with _UPLOAD_LOCK:
        if len(_UPLOAD_CACHE) >= _UPLOAD_MAX:
            _UPLOAD_CACHE.pop(next(iter(_UPLOAD_CACHE)))
        _UPLOAD_CACHE[sub_id] = entry
    return sub_id, entry


_EVO_CACHE: dict[str, dict] = {}


def _get_evolution(pdb_id: str) -> dict | None:
    pid = _norm_id(pdb_id)
    if pid in _EVO_CACHE:
        return _EVO_CACHE[pid]
    _text, structure, _meta = _load_structure(pid)
    uniprots = _get_uniprots(pid)
    evo = evolution.analyze(structure, uniprots)
    if evo is None or evo.get("available") is False:
        _EVO_CACHE[pid] = evo or {
            "available": False,
            "reason": "Conservation analysis is unavailable for this structure.",
        }
        return _EVO_CACHE[pid]
    evo["pocket_conservation"] = evolution.annotate_pockets(evo, _get_pockets(pid))
    evo.pop("_cons_by_key", None)     # internal only
    evo.pop("_hub_keys", None)        # internal only
    evo.pop("_divergent_keys", None)  # internal only
    _EVO_CACHE[pid] = evo
    return evo


_UNIPROT_CACHE: dict[str, list] = {}


def _get_uniprots(pdb_id: str) -> list:
    pid = _norm_id(pdb_id)
    if pid in _UNIPROT_CACHE:
        return _UNIPROT_CACHE[pid]
    try:
        accs = rcsb.fetch_uniprot_accessions(pid)
    except FetchError:
        accs = []
    _UNIPROT_CACHE[pid] = accs
    return accs


def _get_pockets(pdb_id: str) -> list:
    pid = _norm_id(pdb_id)
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


def _get_grid(pdb_id, structure, center):
    """Return a cached docking grid for (pdb_id, site center), building if needed."""
    pid = _norm_id(pdb_id)
    key = (pid, tuple(round(c, 1) for c in center))
    with _GRID_LOCK:
        grid = _GRID_CACHE.get(key)
    if grid is not None:
        return grid
    grid = docking.build_grid(structure, center)
    with _GRID_LOCK:
        if len(_GRID_CACHE) >= _GRID_MAX:
            _GRID_CACHE.pop(next(iter(_GRID_CACHE)))
        _GRID_CACHE[key] = grid
    return grid


def _methods_block(meta, site_label, center, search, ligand=None):
    """Assemble a reproducibility/methods record for a docking or screening run.

    Follows the docking literature's reporting recommendations: tool version,
    receptor identity, box definition, grid spacing, scoring model, search
    settings + random seed, ligand source/flexibility, and prep assumptions.
    """
    block = {
        "tool": f"SnaCleX v{SNACLEX_VERSION}",
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
    benchmark = provenance.docking_benchmark()
    if benchmark is not None:
        block["benchmark"] = benchmark
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


_CHEMCOMP_CACHE: dict[str, dict] = {}


def _resolve_chem_components(codes: list[str]) -> dict:
    """Return {CODE: ccd_info} from the RCSB Chemical Component Dictionary.

    Cached per code; missing/unresolvable codes are cached as ``{}`` so we don't
    refetch them. Any upstream failure degrades to empty (names just won't show).
    """
    want = [c for c in {(c or "").upper() for c in codes} if c and c not in _CHEMCOMP_CACHE]
    if want:
        try:
            fetched = rcsb.fetch_chem_components(want)
        except FetchError:
            fetched = {}
        for code in want:
            _CHEMCOMP_CACHE[code] = fetched.get(code) or {}
    return {c.upper(): _CHEMCOMP_CACHE.get(c.upper()) or {} for c in codes}


def _components_json(structure) -> list[dict]:
    # Resolve real chemical names/chemistry for bound components (ligands, ions,
    # metals) in one batched CCD lookup, so the UI can show "Imatinib" not "STI".
    ccd = _resolve_chem_components([c.res_name for c in structure.components])
    out = []
    for i, c in enumerate(structure.components):
        info = ccd.get((c.res_name or "").upper()) or {}
        out.append(
            {
                "index": i,
                "label": c.label,
                "res_name": c.res_name,
                "chain": c.chain,
                "res_seq": c.res_seq,
                "kind": c.kind,
                "atom_count": len(c.atoms),
                "chem_name": info.get("name"),
                "formula": info.get("formula"),
                "weight": info.get("formula_weight"),
                "smiles": info.get("smiles"),
                "synonyms": info.get("synonyms") or [],
            }
        )
    return out


# ---------------------------------------------------------------------------
# Job runners — the compute behind POST /api/jobs. They take a params dict,
# return a JSON-serializable result, and raise FetchError/ValueError on bad
# input (the JobManager records the message as the job error).
# ---------------------------------------------------------------------------

def _site_params(params):
    """Normalize 'comp'/'pocket' job params to the strings _resolve_dock_site wants."""
    comp = params.get("comp")
    pocket = params.get("pocket")
    return ("" if comp is None else str(comp), "" if pocket is None else str(pocket))


def run_dock_job(params: dict) -> dict:
    pdb_id = params.get("pdb") or ""
    chem = clean_text(params.get("chem") or "")
    comp_raw, pocket_raw = _site_params(params)
    if not pdb_id or not chem or (comp_raw == "" and pocket_raw == ""):
        raise ValueError("Need 'pdb', 'chem', and a site ('comp' or 'pocket')")

    _text, structure, meta = _load_structure(pdb_id)
    center, site_label, ref_component = _resolve_dock_site(
        pdb_id, structure, comp_raw, pocket_raw
    )

    compound = pubchem.lookup_compound(chem)
    if not compound.get("cid"):
        raise ValueError(f"Could not resolve chemical '{chem}'")
    lig = pubchem.fetch_3d_atoms(compound["cid"])
    res_name = (compound.get("molecular_formula") or "LIG")[:3].upper()

    grid = _get_grid(pdb_id, structure, center)
    pose = docking.dock_with_grid(grid, lig["atoms"], center)
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

    return {
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
            meta, site_label, pose["center"], pose["search"],
            ligand={
                "cid": compound.get("cid"),
                "conformer": lig["source"],
                "n_heavy_atoms": len(lig["atoms"]),
            },
        ),
    }


def run_screen_job(params: dict) -> dict:
    import re

    pdb_id = params.get("pdb") or ""
    chems_raw = clean_text(params.get("chems") or "", max_len=MAX_CHEMS_LEN)
    comp_raw, pocket_raw = _site_params(params)
    if not pdb_id or not chems_raw:
        raise ValueError("Need 'pdb' and 'chems' (comma-separated)")

    tokens = [
        t.strip()[:MAX_QUERY_LEN]
        for t in re.split(r"[,;\n]+", chems_raw)
        if t.strip()
    ]
    # De-duplicate, preserve order, cap to keep runtime bounded.
    seen = set()
    chem_list = []
    for t in tokens:
        if t.lower() not in seen:
            seen.add(t.lower())
            chem_list.append(t)
    chem_list = chem_list[:10]
    if not chem_list:
        raise ValueError("No chemicals parsed from 'chems'")

    _text, structure, meta = _load_structure(pdb_id)
    center, site_label, _ref = _resolve_dock_site(
        pdb_id, structure, comp_raw, pocket_raw
    )

    # Build the scoring grid once (cached), dock every ligand against it.
    grid = _get_grid(pdb_id, structure, center)
    results = []
    for token in chem_list:
        try:
            compound = pubchem.lookup_compound(token)
            if not compound.get("cid"):
                raise FetchError("not found in PubChem")
            lig = pubchem.fetch_3d_atoms(compound["cid"])
            pose = docking.dock_with_grid(grid, lig["atoms"], center, seeds=160)
            res_name = (compound.get("molecular_formula") or "LIG")[:3].upper()
            comp = docking.pose_to_component(pose, res_name)
            profile = interactions.profile_component(structure, comp)
            top = [
                f"{r['res_name']}{r['res_seq']}"
                for r in profile["contact_residues"][:3]
            ]
            results.append({
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
            })
        except (FetchError, ValueError) as exc:
            results.append({"query": token, "error": str(exc)})

    # Rank: best (lowest) score first; failures last.
    results.sort(
        key=lambda r: r.get("score", float("inf")) if "error" not in r else float("inf")
    )
    for rank, r in enumerate(results, start=1):
        if "error" not in r:
            r["rank"] = rank

    methods = _methods_block(
        meta, site_label, center, {"seeds": 160, "mc_steps": 40, "random_seed": 0}
    )
    return {
        "pdb_id": _norm_id(pdb_id),
        "site": site_label,
        "count": len(results),
        "results": results,
        "methods": methods,
    }


# Curated known cases for Benchmark Mode (research-credibility checks). The
# server falls back to the largest ligand if the named residue isn't present,
# so a slightly-off resname still benchmarks the main ligand.
BENCHMARK_CASES = [
    {"pdb": "1HSG", "ligand": "MK1", "name": "indinavir",
     "site": "HIV-1 protease active site"},
    {"pdb": "3PTB", "ligand": "BEN", "name": "benzamidine",
     "site": "trypsin S1 pocket"},
    {"pdb": "4DFR", "ligand": "MTX", "name": "methotrexate",
     "site": "dihydrofolate reductase folate site"},
]


def _resolve_benchmark_component(structure, params):
    comp_raw = params.get("comp")
    if comp_raw not in (None, ""):
        idx = int(comp_raw)
        if idx < 0 or idx >= len(structure.components):
            raise ValueError("Component index out of range")
        return structure.components[idx]
    ligs = structure.ligand_components
    if not ligs:
        raise ValueError("No bound ligand in this structure to benchmark.")
    wanted = clean_text(params.get("ligand") or "", max_len=8).upper()
    if wanted:
        matches = [c for c in structure.components if c.res_name.upper() == wanted]
        if matches:
            return max(matches, key=lambda c: len(c.atoms))
    # Fall back to the largest ligand component.
    return max(ligs, key=lambda c: len(c.atoms))


def run_benchmark_job(params: dict) -> dict:
    pdb_id = params.get("pdb") or ""
    if not pdb_id:
        raise ValueError("Need 'pdb'")
    _text, structure, meta = _load_structure(pdb_id)
    comp = _resolve_benchmark_component(structure, params)
    result = benchmark.benchmark_case(structure, comp)
    result["receptor"] = {"pdb_id": meta.get("pdb_id"), "title": meta.get("title")}
    result["methods"] = _methods_block(
        meta, comp.label, docking.component_center(comp), result["search"]
    )
    return result


def _default_interface_groups(structure):
    """Heuristic chain grouping: the two smallest chains vs the rest.

    A reasonable zero-config default (e.g. antibody H+L, or a short peptide,
    against a larger antigen); the user can always override with explicit a/b.
    """
    counts: dict[str, int] = {c: 0 for c in structure.chains}
    for a in structure.protein_atoms:
        counts[a.chain] = counts.get(a.chain, 0) + 1
    chains = sorted(counts, key=lambda c: counts[c])
    if len(chains) < 2:
        return [], []
    if len(chains) == 2:
        return [chains[0]], [chains[1]]
    return chains[:2], chains[2:]


_GENOMIC_RE = re.compile(
    r"^(?:chr)?([0-9]{1,2}|[XYM]|MT)[:\-]([0-9]+)[:\-]([ACGT]+)[>\-/]([ACGT]+)$",
    re.IGNORECASE,
)


def _residue_pocket_index(pockets_list):
    idx = {}
    for p in pockets_list:
        for r in p.get("lining_residues", []):
            idx.setdefault((r["chain"], r["res_seq"]), p["index"])
    return idx


def _residue_conservation(evo):
    out = {}
    if evo and evo.get("residues"):
        for r in evo["residues"]:
            if r.get("conservation") is not None:
                out[(r["chain"], r["res_seq"])] = r["conservation"]
    return out


def run_variants_job(params: dict) -> dict:
    """Map a list of variants onto a structure and annotate them.

    Protein-position input (R273H / p.Arg273His / position) is mapped to
    residues now; genomic input (chr:pos ref>alt) gets a best-effort TOPMed/BRAVO
    allele frequency. When the structure is HLA, mapped variants are also
    classified against the peptide-binding groove. Optional ESM variant-effect
    scores are attached alongside when a Forge token is configured.
    """
    pdb_id = params.get("pdb") or ""
    raw_variants = params.get("variants") or []
    acc_override = clean_text(str(params.get("uniprot") or ""), max_len=20)
    if not pdb_id:
        raise ValueError("Need 'pdb'")
    if not isinstance(raw_variants, list) or not raw_variants:
        raise ValueError("Need a non-empty 'variants' list")

    variant_strs = []
    for v in raw_variants[:200]:
        s = clean_text(str(v), max_len=32)
        if s:
            variant_strs.append(s)

    # Genomic-format inputs go to the BRAVO frequency path; the rest are treated
    # as protein-position variants for structural mapping.
    protein_inputs, genomic = [], []
    for s in variant_strs:
        if _GENOMIC_RE.match(s):
            genomic.append(s)
        else:
            protein_inputs.append(s)

    _text, structure, meta = _load_structure(pdb_id)
    accs = [acc_override] if acc_override else _get_uniprots(pdb_id)

    uniprot_seq, acc_used = "", None
    for acc in accs:
        try:
            seq = variants.fetch_uniprot_sequence(acc)
        except FetchError:
            seq = ""
        if seq:
            uniprot_seq, acc_used = seq, acc
            break

    mapping = variants.annotate(structure, protein_inputs, uniprot_seq)

    # Genomic input -> protein consequence (Ensembl VEP, env-gated) -> residue.
    # Merged into the same variant list so it gets the same enrichment + HLA
    # classification below. When VEP is off, genomic variants keep only their
    # BRAVO frequency (added later) without a structural location.
    if genomic and vep.available():
        seq_cache = {acc_used: uniprot_seq} if acc_used else {}
        for s in genomic:
            m = _GENOMIC_RE.match(s)
            cons = vep.annotate_one(m.group(1), m.group(2), m.group(3), m.group(4))
            entry = {"input": s, "source": "genomic"}
            if not cons:
                entry.update(mapped=False, reason="no protein-coding consequence (VEP)")
                mapping["variants"].append(entry)
                continue
            entry.update(gene=cons["gene"], consequence=cons["consequence"],
                         wt=cons["wt_aa"], mut=cons["mut_aa"],
                         position=cons["protein_position"], uniprot=cons["uniprot"])
            if cons["uniprot"] not in accs:
                entry.update(mapped=False,
                             reason=f"affects {cons['gene'] or cons['uniprot']}, "
                                    "not this structure's protein")
                mapping["variants"].append(entry)
                continue
            acc = cons["uniprot"]
            if acc not in seq_cache:
                try:
                    seq_cache[acc] = variants.fetch_uniprot_sequence(acc)
                except FetchError:
                    seq_cache[acc] = ""
            sub = f"{cons['wt_aa']}{cons['protein_position']}{cons['mut_aa']}"
            mapped = variants.annotate(structure, [sub], seq_cache[acc])["variants"][0]
            mapped.update(input=s, source="genomic", gene=cons["gene"],
                          consequence=cons["consequence"])
            mapping["variants"].append(mapped)
        mapping["input_count"] = len(mapping["variants"])
        mapping["mapped_count"] = sum(1 for r in mapping["variants"] if r.get("mapped"))

    # Structural context for mapped residues (pocket membership + conservation).
    pocket_idx = _residue_pocket_index(_get_pockets(pdb_id))
    try:
        cons = _residue_conservation(_get_evolution(pdb_id))
    except FetchError:
        cons = {}
    for v in mapping["variants"]:
        if v.get("mapped"):
            key = (v["chain"], v["res_seq"])
            v["pocket"] = pocket_idx.get(key)
            v["conservation"] = cons.get(key)

    # HLA groove classification (ties variants to the peptide-binding groove).
    hla_block = None
    detection = hla.detect(structure, title=meta.get("title"), uniprots=accs)
    if detection.get("is_hla"):
        groove = hla.analyze_groove(structure, detection)
        hla_block = {
            "detection": detection,
            "groove": groove,
            "variant_classes": hla.classify_variants(groove, mapping["variants"]),
        }

    # TOPMed/BRAVO allele frequencies for genomic-format inputs (best-effort).
    population = []
    for s in genomic:
        m = _GENOMIC_RE.match(s)
        chrom, pos, ref, alt = m.group(1), m.group(2), m.group(3), m.group(4)
        freq = bravo.variant_frequency(chrom, pos, ref, alt)
        if freq.get("available"):
            freq["rarity"] = bravo.classify_frequency(freq.get("allele_freq"))
        freq["input"] = s
        population.append(freq)

    # ESM variant-effect scores alongside (only if a Forge token is set).
    esm = (models_esm.score_variants(uniprot_seq, protein_inputs)
           if models_esm.available() and uniprot_seq
           else {"available": False, "reason": "ESM scoring not enabled (no ESM_API_KEY)"})

    return {
        "pdb_id": meta.get("pdb_id") or pdb_id,
        "uniprot": acc_used,
        "uniprot_candidates": accs,
        "mapping": mapping,
        "hla": hla_block,
        "population": population,
        "esm": esm,
        "vep_enabled": vep.available(),
        "methods": provenance.variant_methods(),
    }


_JOB_RUNNERS = {
    "dock": run_dock_job,
    "screen": run_screen_job,
    "benchmark": run_benchmark_job,
    "variants": run_variants_job,
}


class Handler(BaseHTTPRequestHandler):
    server_version = "SnaCleX/0.1"

    def log_message(self, fmt, *args):  # quieter console
        pass

    # ---- helpers ------------------------------------------------------
    def _common_headers(self):
        """Security headers emitted on every response.

        No Access-Control-Allow-Origin is sent: the API is deliberately
        same-origin only (the browser talks to /api/*, which proxies the
        upstream scientific services server-side).
        """
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("X-Frame-Options", "DENY")
        self.send_header("Content-Security-Policy", _CSP)
        if self.headers.get("X-Forwarded-Proto", "").lower() == "https":
            self.send_header(
                "Strict-Transport-Security", "max-age=31536000; includeSubDomains"
            )

    def _client_ip(self):
        fwd = self.headers.get("X-Forwarded-For", "")
        if fwd:
            return fwd.split(",")[0].strip()
        return self.client_address[0] if self.client_address else "unknown"

    def _send_json(self, payload, status=200, retry_after=None):
        self._status = status
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        if retry_after is not None:
            self.send_header("Retry-After", str(int(retry_after) + 1))
        self._common_headers()
        self.end_headers()
        self.wfile.write(body)

    def _send_error_json(self, message, status=400, retry_after=None):
        self._send_json({"error": message}, status=status, retry_after=retry_after)

    def _serve_static(self, path):
        if path in ("/", ""):
            path = "/index.html"
        rel = path.lstrip("/")
        full = os.path.normpath(os.path.join(WEB_DIR, rel))
        if not full.startswith(WEB_DIR) or not os.path.isfile(full):
            self._status = 404
            self.send_error(404, "Not found")
            return
        self._status = 200
        ext = os.path.splitext(full)[1].lower()
        ctype = _CONTENT_TYPES.get(ext, "application/octet-stream")
        with open(full, "rb") as fh:
            body = fh.read()
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self._common_headers()
        self.end_headers()
        self.wfile.write(body)

    # ---- access logging ----------------------------------------------
    def _access_log(self, method, path):
        dur_ms = (time.monotonic() - getattr(self, "_t0", time.monotonic())) * 1000
        status = getattr(self, "_status", None) or "-"
        # Job-status polling is high-frequency and low-signal — keep it at DEBUG.
        level = logging.DEBUG if path.startswith("/api/jobs/") else logging.INFO
        log.log(level, "%s %s -> %s %.0fms ip=%s",
                method, path, status, dur_ms, _ip_hash(self._client_ip()))

    # ---- routing ------------------------------------------------------
    def do_GET(self):
        self._t0 = time.monotonic()
        self._status = None
        path = urlparse(self.path).path
        try:
            parsed = urlparse(self.path)
            self._dispatch_get(parsed.path, parse_qs(parsed.query))
        finally:
            self._access_log("GET", path)

    def _dispatch_get(self, path, qs):
        # Job-status polling is cheap and frequent — exempt it from rate limits
        # (job ids are unguessable, so this isn't an abuse vector).
        if path.startswith("/api/jobs/"):
            return self._guarded_route(path, qs)

        if path.startswith("/api/"):
            ip = self._client_ip()
            ok, retry = _IP_LIMITER.allow(ip)
            if not ok:
                return self._send_error_json(
                    "Too many requests; please slow down.",
                    status=429, retry_after=retry,
                )
            if path in EXPENSIVE_ENDPOINTS:
                ok, retry = _HEAVY_LIMITER.allow(ip)
                if not ok:
                    return self._send_error_json(
                        "This analysis is rate-limited; try again shortly.",
                        status=429, retry_after=retry,
                    )
                if not _HEAVY_SEM.acquire(blocking=False):
                    return self._send_error_json(
                        "Server is busy running analyses; please retry in a moment.",
                        status=503, retry_after=5,
                    )
                try:
                    return self._guarded_route(path, qs)
                finally:
                    _HEAVY_SEM.release()

        return self._guarded_route(path, qs)

    def _guarded_route(self, path, qs):
        try:
            return self._route(path, qs)
        except FetchError as exc:
            return self._send_error_json(str(exc), status=502)
        except Exception:  # noqa: BLE001
            # Log the traceback server-side; return a generic message to the
            # client so internal details aren't leaked.
            log.exception("Unhandled error serving %s", path)
            return self._send_error_json(
                "Internal error — please retry.", status=500
            )

    def _route(self, path, qs):
        if path == "/api/analyze":
            return self._api_analyze(qs)
        if path == "/api/interactions":
            return self._api_interactions(qs)
        if path == "/api/interface":
            return self._api_interface(qs)
        if path == "/api/nucleic":
            return self._api_nucleic(qs)
        if path == "/api/hla":
            return self._api_hla(qs)
        if path == "/api/hla/cases":
            return self._send_json({"cases": hla.curated_cases()})
        if path == "/api/expression":
            return self._api_expression(qs)
        if path == "/api/esm":
            return self._send_json(models_esm.config())
        if path == "/api/chemical":
            return self._api_chemical(qs)
        if path == "/api/pockets":
            return self._api_pockets(qs)
        if path == "/api/evolution":
            return self._api_evolution(qs)
        if path == "/api/search":
            return self._api_search(qs)
        if path == "/api/version":
            return self._api_version(qs)
        if path == "/api/docs":
            return self._send_json(apidocs.contract())
        if path == "/api/benchmark/cases":
            return self._send_json({"cases": BENCHMARK_CASES})
        if path.startswith("/api/jobs/"):
            return self._api_job_status(path)
        return self._serve_static(path)

    # ---- POST endpoints ----------------------------------------------
    def do_POST(self):
        self._t0 = time.monotonic()
        self._status = None
        path = urlparse(self.path).path
        try:
            self._dispatch_post(path)
        except FetchError as exc:
            self._send_error_json(str(exc), status=502)
        except Exception:  # noqa: BLE001
            log.exception("Unhandled error serving POST %s", path)
            self._send_error_json("Internal error — please retry.", status=500)
        finally:
            self._access_log("POST", path)

    def _dispatch_post(self, path):
        if path not in ("/api/jobs", "/api/upload"):
            return self._send_error_json("Not found", status=404)

        ip = self._client_ip()
        ok, retry = _IP_LIMITER.allow(ip)
        if not ok:
            return self._send_error_json(
                "Too many requests; please slow down.", status=429, retry_after=retry
            )
        ok, retry = _HEAVY_LIMITER.allow(ip)
        if not ok:
            return self._send_error_json(
                "This analysis is rate-limited; try again shortly.",
                status=429, retry_after=retry,
            )

        if path == "/api/upload":
            return self._guarded_upload()

        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0 or length > 8192:
            return self._send_error_json("Missing or oversized request body")
        try:
            payload = json.loads(self.rfile.read(length))
        except (ValueError, OSError):
            return self._send_error_json("Invalid JSON body")
        if not isinstance(payload, dict):
            return self._send_error_json("Body must be a JSON object")

        kind = payload.get("kind")
        if kind not in _JOB_RUNNERS:
            return self._send_error_json(f"Unknown job kind '{kind}'")
        params = payload.get("params") or {}
        if not isinstance(params, dict):
            return self._send_error_json("'params' must be an object")

        job_id = JOBS.submit(_JOB_RUNNERS[kind], params)
        return self._send_json({"job_id": job_id, "status": "queued"}, status=202)

    def _guarded_upload(self):
        try:
            return self._api_upload()
        except FetchError as exc:  # size guard et al. — surface the message
            return self._send_error_json(str(exc), status=413)
        except Exception as exc:  # noqa: BLE001
            return self._send_error_json(f"Could not parse structure: {exc}")

    def _api_upload(self):
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0:
            return self._send_error_json("Empty upload")
        if length > MAX_UPLOAD_BYTES:
            return self._send_error_json(
                f"Structure too large (max {MAX_UPLOAD_BYTES // 1_000_000} MB)",
                status=413,
            )
        raw = self.rfile.read(length).decode("utf-8", errors="replace")

        structure = pdbparse.parse_structure(raw)  # PDB or mmCIF, auto-detected
        if len(structure.protein_atoms) < 10:
            return self._send_error_json(
                "No protein atoms found — is this a valid PDB/mmCIF file?"
            )
        _check_structure_size(structure, "This file")

        # The viewer reads PDB format, so re-serialize mmCIF uploads.
        viewer_text = pdbparse.to_pdb(structure) if "_atom_site." in raw else raw

        upload_id = "UL" + uuid.uuid4().hex[:10]
        meta = {"pdb_id": upload_id, "title": "Uploaded structure", "uploaded": True}
        with _UPLOAD_LOCK:
            if len(_UPLOAD_CACHE) >= _UPLOAD_MAX:
                _UPLOAD_CACHE.pop(next(iter(_UPLOAD_CACHE)))
            _UPLOAD_CACHE[upload_id] = (viewer_text, structure, meta)

        return self._send_json({
            "upload_id": upload_id,
            "metadata": meta,
            "chains": structure.chains,
            "nucleic_chains": structure.nucleic_chains,
            "protein_atom_count": len(structure.protein_atoms),
            "nucleic_atom_count": len(structure.nucleic_atoms),
            "components": _components_json(structure),
            "pdb_data": viewer_text,
        })

    def _api_job_status(self, path):
        job_id = path[len("/api/jobs/"):]
        job = JOBS.status(job_id)
        if job is None:
            return self._send_error_json("Unknown or expired job", status=404)
        out = {"job_id": job_id, "status": job["status"]}
        if job["status"] == "done":
            out["result"] = job["result"]
        elif job["status"] == "error":
            out["error"] = job["error"]
        return self._send_json(out)

    # ---- API endpoints ------------------------------------------------
    def _analyze_payload(self, sid, structure, meta, text):
        return self._send_json({
            "id": sid,
            "metadata": meta,
            "chains": structure.chains,
            "nucleic_chains": structure.nucleic_chains,
            "protein_atom_count": len(structure.protein_atoms),
            "nucleic_atom_count": len(structure.nucleic_atoms),
            "components": _components_json(structure),
            "pdb_data": text,
        })

    def _api_analyze(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        chain = clean_text((qs.get("chain") or [""])[0], max_len=8)
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")

        # Uploaded structures / per-chain subsets are already small and cached.
        if pdb_id in _UPLOAD_CACHE:
            text, structure, meta = _load_structure(pdb_id)
            return self._analyze_payload(pdb_id, structure, meta, text)

        raw, structure, meta = _load_full(pdb_id)
        pid = meta.get("pdb_id") or rcsb.normalize_pdb_id(pdb_id)

        if chain:
            if chain not in structure.chains:
                return self._send_error_json(
                    f"Chain '{chain}' not found in {pid}.", status=404
                )
            sub_id, (text, sub, sub_meta) = _cache_subset(pid, chain, structure, meta)
            return self._analyze_payload(sub_id, sub, sub_meta, text)

        if len(structure.atoms) > MAX_STRUCTURE_ATOMS:
            counts = {}
            for a in structure.atoms:
                counts[a.chain] = counts.get(a.chain, 0) + 1
            chains = [{"chain": c, "atom_count": counts.get(c, 0)}
                      for c in structure.chains]
            return self._send_json({
                "too_large": True,
                "pdb_id": pid,
                "n_atoms": len(structure.atoms),
                "limit": MAX_STRUCTURE_ATOMS,
                "chains": chains,
            })

        text = pdbparse.to_pdb(structure) if "_atom_site." in raw else raw
        return self._analyze_payload(pid, structure, meta, text)

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

    def _api_interface(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        a = clean_text((qs.get("a") or [""])[0], max_len=64)
        b = clean_text((qs.get("b") or [""])[0], max_len=64)
        heavy_ch = clean_text((qs.get("heavy") or [""])[0], max_len=8) or None
        light_ch = clean_text((qs.get("light") or [""])[0], max_len=8) or None
        scheme = clean_text((qs.get("scheme") or [""])[0], max_len=16) or "kabat"
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        _text, structure, meta = _load_structure(pdb_id)
        chains_a = [c for c in a.split(",") if c] if a else []
        chains_b = [c for c in b.split(",") if c] if b else []
        if not chains_a or not chains_b:
            # Default heuristic: the two shortest chains (e.g. antibody H+L) vs
            # the rest, so the user gets a sensible interface without picking.
            chains_a, chains_b = _default_interface_groups(structure)
        if not chains_a or not chains_b:
            return self._send_error_json(
                "Need at least two protein chains to define an interface.", status=404
            )
        profile = interactions.profile_interface(structure, chains_a, chains_b)

        # CDR annotation — explicit heavy/light takes priority; fall back to
        # auto-detection when neither is provided but chains_a looks antibody-like.
        if heavy_ch or light_ch:
            profile = antibody.annotate_paratope(profile, heavy_ch, light_ch, scheme)
        else:
            det = antibody.detect_vhvl_chains(structure)
            if det["confidence"] in ("high", "medium") and (det["heavy"] or det["light"]):
                # Only annotate if the detected chains are actually in chains_a.
                h = det["heavy"] if det["heavy"] in set(chains_a) else None
                l = det["light"] if det["light"] in set(chains_a) else None
                if h or l:
                    profile = antibody.annotate_paratope(profile, h, l, scheme)
                    profile["antibody_auto_detected"] = True

        return self._send_json({
            "profile": profile,
            "metadata": meta,
            "methods": provenance.interface_methods(),
        })

    def _api_nucleic(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        chain = clean_text((qs.get("chain") or [""])[0], max_len=8) or None
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        _text, structure, meta = _load_structure(pdb_id)
        if not structure.nucleic_atoms:
            return self._send_json({
                "profile": None,
                "metadata": meta,
                "available": False,
                "reason": "No nucleic-acid (DNA/RNA) chains in this structure.",
            })
        profile = interactions.profile_nucleic_interface(structure, chain)
        return self._send_json({
            "profile": profile,
            "metadata": meta,
            "available": True,
            "methods": provenance.interface_methods(),
        })

    def _api_hla(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        _text, structure, meta = _load_structure(pdb_id)
        uniprots = _get_uniprots(pdb_id)
        detection = hla.detect(structure, title=meta.get("title"), uniprots=uniprots)
        groove = hla.analyze_groove(structure, detection) if detection["is_hla"] else None
        return self._send_json({
            "metadata": meta,
            "detection": detection,
            "groove": groove,
            "methods": provenance.hla_methods(),
        })

    def _api_expression(self, qs):
        gene = clean_text((qs.get("gene") or [""])[0], max_len=32)
        if not gene:
            return self._send_error_json("Missing 'gene' parameter")
        result = cellxgene.gene_expression(gene)
        result["note"] = (
            "Bulk single-cell expression context from CZ CELL×GENE — which cell "
            "types express this gene across human tissues. Research-only; not "
            "patient-specific. The summary requires SNACLEX_ENABLE_CELLXGENE; the "
            "deep link is always available."
        )
        return self._send_json(result)

    def _api_chemical(self, qs):
        query = clean_text((qs.get("q") or [""])[0])
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

    def _api_pockets(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        found = _get_pockets(pdb_id)
        return self._send_json({
            "pockets": found,
            "count": len(found),
            "methods": provenance.pocket_methods(),
        })

    def _api_evolution(self, qs):
        pdb_id = (qs.get("pdb") or [""])[0]
        if not pdb_id:
            return self._send_error_json("Missing 'pdb' parameter")
        evo = _get_evolution(pdb_id)
        if evo is None or evo.get("available") is False:
            return self._send_json(
                evo or {"available": False, "reason": "Conservation analysis unavailable."}
            )
        return self._send_json({
            "available": True,
            **evo,
            "methods": provenance.evolution_methods(),
        })

    def _api_search(self, qs):
        query = clean_text((qs.get("q") or [""])[0])
        if not query:
            return self._send_error_json("Missing 'q' parameter")
        return self._send_json({"results": rcsb.search_by_name(query, limit=10)})

    def _api_version(self, qs):
        return self._send_json({
            "name": "SnaCleX",
            "version": SNACLEX_VERSION,
            "research_only": True,
        })


def main():
    # Cloud hosts (Render/Railway/Fly/etc.) inject the port via $PORT and need
    # the server bound to all interfaces. Locally these default to 8010/127.0.0.1.
    env_port = os.environ.get("PORT")
    default_port = int(env_port) if env_port and env_port.isdigit() else 8010
    default_host = "0.0.0.0" if env_port else "127.0.0.1"

    parser = argparse.ArgumentParser(description="SnaCleX web server")
    parser.add_argument("--port", type=int, default=default_port)
    parser.add_argument("--host", default=default_host)
    args = parser.parse_args()

    configure_logging()
    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    log.info("SnaCleX v%s listening on %s:%s", SNACLEX_VERSION, args.host, args.port)
    print(f"SnaCleX running at http://{args.host}:{args.port}")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        httpd.shutdown()


if __name__ == "__main__":
    main()
