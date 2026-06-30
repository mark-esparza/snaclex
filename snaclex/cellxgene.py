"""CELL×GENE expression context (supplemental, best-effort).

Reports which tissues / immune cell types express a gene, to contextualize an
HLA or antigen finding, from CZ CELL×GENE Discover. Two tiers:

  * ``discover_link(gene)`` — a deep link into the CELL×GENE Gene Expression app.
    Always available, no network; this is the reliable baseline.
  * ``gene_expression(gene)`` — a best-effort summary of the top expressing cell
    types via the public "Where's My Gene" (WMG) API. OFF unless
    ``SNACLEX_ENABLE_CELLXGENE`` is set (external dependency; the remote env's
    network policy may block it). Resolves the gene symbol to an Ensembl gene id
    (Ensembl REST) first. Behind an injectable seam; any failure degrades to the
    link only.

Research-only: bulk single-cell expression context, not patient-specific data.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

from .http_util import FetchError, fetch_json

_ENV_FLAG = "SNACLEX_ENABLE_CELLXGENE"
_WMG_BASE = "https://api.cellxgene.cziscience.com/wmg/v2"
_HUMAN = "NCBITaxon:9606"


def available() -> bool:
    v = (os.environ.get(_ENV_FLAG) or "").strip().lower()
    return v not in ("", "0", "off", "false", "no")


def _default_fetcher(url: str, body=None):
    """GET (cached) when ``body`` is None, else POST JSON. Raises FetchError."""
    if body is None:
        return fetch_json(url)
    data = json.dumps(body).encode("utf-8")
    req = urllib.request.Request(
        url, data=data,
        headers={"Content-Type": "application/json",
                 "User-Agent": "SnaCleX/0.1 (research tool)"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8", errors="replace"))
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        raise FetchError(f"CELL×GENE POST failed: {exc}") from exc


def discover_link(gene: str) -> str:
    """Deep link into the CELL×GENE Gene Expression app for a gene symbol."""
    from urllib.parse import quote
    return f"https://cellxgene.cziscience.com/gene-expression?genes={quote(gene or '')}"


def _resolve_ensembl_id(symbol: str, fetcher) -> str | None:
    """Gene symbol -> Ensembl gene id (ENSG...) via Ensembl xrefs."""
    url = (f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
           "?content-type=application/json")
    try:
        rows = fetcher(url)
    except FetchError:
        return None
    for r in rows or []:
        gid = r.get("id")
        if gid and str(gid).startswith("ENSG"):
            return gid
    return None


def _top_cell_types(payload, ensembl_id, limit=12) -> list:
    """Parse a WMG query response into ranked {cell_type, mean_expr, pct_cells}.

    Tolerant of schema drift: returns [] if the expected shape isn't present.
    """
    summary = (payload or {}).get("expression_summary") or {}
    gene_block = summary.get(ensembl_id) or {}
    labels = ((payload or {}).get("term_id_labels") or {}).get("cell_types") or {}
    name_of = {}
    if isinstance(labels, dict):
        name_of = {k: (v.get("name") if isinstance(v, dict) else v) for k, v in labels.items()}
    elif isinstance(labels, list):
        for it in labels:
            if isinstance(it, dict) and it.get("cell_type_ontology_term_id"):
                name_of[it["cell_type_ontology_term_id"]] = it.get("name")

    rows = []
    for _tissue, cells in (gene_block.items() if isinstance(gene_block, dict) else []):
        if not isinstance(cells, dict):
            continue
        for cell_id, stats in cells.items():
            if not isinstance(stats, dict):
                continue
            me = stats.get("me")  # mean expression
            if me is None:
                continue
            rows.append({
                "cell_type": name_of.get(cell_id, cell_id),
                "mean_expr": round(float(me), 3),
                "pct_cells": stats.get("pc"),
            })
    rows.sort(key=lambda r: -(r["mean_expr"] or 0))
    return rows[:limit]


def gene_expression(gene: str, fetcher=None) -> dict:
    """Best-effort top-expressing cell types for a gene symbol.

    Always returns the deep link; ``cell_types`` is populated only when enabled,
    the gene resolves, and the WMG API is reachable.
    """
    gene = (gene or "").strip()
    out = {"gene": gene, "link": discover_link(gene), "available": False}
    if not gene:
        out["reason"] = "no gene symbol"
        return out
    if fetcher is None and not available():
        out["reason"] = "CELL×GENE summary not enabled (set SNACLEX_ENABLE_CELLXGENE)"
        return out
    fetcher = fetcher or _default_fetcher

    ensembl_id = _resolve_ensembl_id(gene, fetcher)
    if not ensembl_id:
        out["reason"] = f"could not resolve an Ensembl gene id for '{gene}'"
        return out

    body = {
        "filter": {
            "gene_ontology_term_ids": [ensembl_id],
            "organism_ontology_term_id": _HUMAN,
            "tissue_ontology_term_ids": [],
        },
        "is_rollup": True,
    }
    try:
        payload = fetcher(f"{_WMG_BASE}/query", body)
    except FetchError as exc:
        out["reason"] = f"CELL×GENE WMG unavailable ({exc})"
        return out

    cells = _top_cell_types(payload, ensembl_id)
    if not cells:
        out["reason"] = "no expression summary returned"
        return out
    out.update(available=True, ensembl_id=ensembl_id, cell_types=cells,
               source="CZ CELL×GENE (WMG)")
    return out
