"""Genome variant bridge: overlay ClinVar / gnomAD missense variants on a structure.

This is a genotype-to-structure BRIDGE, not a genome browser: no chromosome
ideograms, no coordinate system, no reads. It resolves the loaded structure to a
UniProt accession, pulls missense variants for that protein, and maps each
variant's UniProt residue position onto the loaded structure's residue numbering
so it can be colored on the 3D model (reusing the Evolution coloring path).

Data source: the EMBL-EBI Proteins API "variation" endpoint, which aggregates
**ClinVar** clinical significance and **gnomAD** population allele frequencies per
UniProt position (and ships the UniProt sequence needed to align onto the
structure). The optional gene/locus card uses NCBI E-utilities. Research-only.
"""

from __future__ import annotations

import datetime

from .evolution import AA3TO1, _nw_align
from .http_util import FetchError, fetch_json

# Bound payloads for pathological proteins with huge variant sets.
MAX_VARIANTS_PER_RESIDUE = 12
MAX_MAPPED_RESIDUES = 4000

# ClinVar significance -> (ordinal severity, canonical label). Higher = worse.
_SIG_RANK = {
    "pathogenic": (5, "pathogenic"),
    "likely pathogenic": (4, "likely pathogenic"),
    "pathogenic/likely pathogenic": (5, "pathogenic"),
    "risk factor": (3, "risk factor"),
    "uncertain significance": (2, "uncertain"),
    "conflicting interpretations of pathogenicity": (2, "conflicting"),
    "not provided": (1, "not provided"),
    "likely benign": (1, "likely benign"),
    "benign/likely benign": (1, "benign"),
    "benign": (0, "benign"),
}


def _unavailable(reason: str) -> dict:
    return {"available": False, "reason": reason}


def _canonical_sig(significances) -> tuple:
    """Reduce a variant's ClinVar significances to (rank, label)."""
    best = (-1, None)
    for s in significances or []:
        t = (s.get("type") or "").strip().lower()
        rank, label = _SIG_RANK.get(t, (2, t or "uncertain"))
        if rank > best[0]:
            best = (rank, label)
    if best[1] is None:
        return (-1, None)
    return best


def _gnomad_af(freqs):
    """Max gnomAD allele frequency across the reported populations, or None."""
    best = None
    for f in freqs or []:
        src = (f.get("source") or "").lower()
        if "gnomad" not in src:
            continue
        val = f.get("frequency")
        try:
            val = float(val)
        except (TypeError, ValueError):
            continue
        if best is None or val > best:
            best = val
    return best


def _parse_variants(features):
    """Extract single-residue missense variants from EBI variation features."""
    out = []
    for f in features or []:
        if f.get("type") != "VARIANT":
            continue
        wt = (f.get("wildType") or "").strip()
        alt = (f.get("alternativeSequence") or f.get("mutatedType") or "").strip()
        if len(wt) != 1 or len(alt) != 1 or not wt.isalpha() or not alt.isalpha():
            continue
        if alt == wt or alt == "*":
            continue
        try:
            pos = int(f.get("begin"))
            if int(f.get("end")) != pos:
                continue
        except (TypeError, ValueError):
            continue

        rank, sig = _canonical_sig(f.get("clinicalSignificances"))
        af = _gnomad_af(f.get("populationFrequencies"))
        if sig is None and af is None:
            continue  # neither ClinVar nor gnomAD says anything — skip noise

        xrefs = {x.get("name"): x.get("id") for x in (f.get("xrefs") or [])}
        out.append({
            "position": pos,
            "ref": wt,
            "alt": alt,
            "clinical_significance": sig,
            "sig_rank": rank,
            "allele_frequency": af,
            "dbsnp": xrefs.get("dbSNP"),
            "clinvar": xrefs.get("ClinVar"),
        })
    return out


def fetch_variation(accession: str) -> dict | None:
    """Fetch ClinVar/gnomAD missense variants + UniProt sequence for an accession."""
    url = f"https://www.ebi.ac.uk/proteins/api/variation/{accession}"
    try:
        data = fetch_json(url)
    except FetchError:
        return None
    if isinstance(data, list):
        data = data[0] if data else {}
    seq = data.get("sequence")
    if not seq:
        return None
    return {
        "accession": data.get("accession") or accession,
        "gene": (data.get("geneName") or "").strip() or None,
        "sequence": seq,
        "variants": _parse_variants(data.get("features")),
    }


def _map_positions(structure, uniprot_seq):
    """Map UniProt 1-based residue positions onto structure residues.

    Aligns each protein chain to the UniProt sequence (reusing the Evolution
    NW aligner) so homomers land variants on every copy. Returns
    {uniprot_pos: [(chain, res_seq, res_name), …]} and the covered residue count.
    """
    by_chain: dict[str, list] = {}
    for a in structure.protein_atoms:
        d = by_chain.setdefault(a.chain, {})
        if a.res_seq not in d:
            d[a.res_seq] = a.res_name
    pos_map: dict[int, list] = {}
    mapped_residues = 0
    for chain, residues in by_chain.items():
        items = sorted(residues.items())  # (res_seq, res_name)
        chain_seq = "".join(AA3TO1.get(rn, "X") for _, rn in items)
        if len(chain_seq) < 15:
            continue
        pairs = _nw_align(chain_seq, uniprot_seq)
        matches = 0
        local = []
        for ti, uj in pairs:
            if ti is None or uj is None:
                continue
            res_seq, res_name = items[ti]
            if chain_seq[ti] == uniprot_seq[uj]:
                matches += 1
            local.append((uj + 1, chain, res_seq, res_name))  # UniProt is 1-based
        # Only trust a chain that clearly is this UniProt protein.
        if not local or matches / len(local) < 0.6:
            continue
        for up_pos, ch, rs, rn in local:
            pos_map.setdefault(up_pos, []).append((ch, rs, rn))
            mapped_residues += 1
    return pos_map, mapped_residues


def analyze(structure, uniprots) -> dict:
    """Full variant-bridge analysis for a loaded structure."""
    if not uniprots:
        return _unavailable(
            "This structure has no UniProt mapping (common for designed, "
            "synthetic or peptide-only entries), so no ClinVar/gnomAD variants "
            "can be cross-referenced."
        )

    data = None
    acc_used = None
    for acc in uniprots:
        data = fetch_variation(acc)
        if data:
            acc_used = acc
            break
    if not data:
        return _unavailable(
            f"No variant record was returned for this protein (UniProt "
            f"{', '.join(uniprots)}) from the EBI Proteins API — the service may "
            f"be unavailable, or the protein has no annotated missense variants."
        )

    variants = data["variants"]
    if not variants:
        return _unavailable(
            f"UniProt {acc_used} has no ClinVar/gnomAD missense variants to map."
        )

    up_seq = data["sequence"]
    pos_map, _covered = _map_positions(structure, up_seq)
    if not pos_map:
        return _unavailable(
            f"Could not align the structure's sequence to UniProt {acc_used}, so "
            f"variant positions cannot be placed on this structure."
        )

    # Group variants onto structure residues.
    residues: dict[tuple, dict] = {}
    mapped_count = 0
    unmapped = 0
    for v in variants:
        targets = pos_map.get(v["position"])
        if not targets:
            unmapped += 1
            continue
        mapped_count += 1
        for chain, res_seq, res_name in targets:
            key = (chain, res_seq)
            entry = residues.setdefault(key, {
                "chain": chain, "res_seq": res_seq, "res_name": res_name,
                "uniprot_pos": v["position"], "variants": [],
                "worst_rank": -1, "pathogenicity": None, "max_af": None,
            })
            if len(entry["variants"]) < MAX_VARIANTS_PER_RESIDUE:
                entry["variants"].append({
                    "ref": v["ref"], "alt": v["alt"],
                    "clinical_significance": v["clinical_significance"],
                    "allele_frequency": v["allele_frequency"],
                    "dbsnp": v["dbsnp"], "clinvar": v["clinvar"],
                })
            if v["sig_rank"] > entry["worst_rank"]:
                entry["worst_rank"] = v["sig_rank"]
                entry["pathogenicity"] = v["clinical_significance"]
            af = v["allele_frequency"]
            if af is not None and (entry["max_af"] is None or af > entry["max_af"]):
                entry["max_af"] = af

    residue_list = sorted(
        residues.values(), key=lambda e: (-e["worst_rank"], -(e["max_af"] or 0))
    )[:MAX_MAPPED_RESIDUES]

    # Pathogenicity summary counts (over mapped residues).
    summary: dict[str, int] = {}
    for e in residue_list:
        label = e["pathogenicity"] or "frequency-only"
        summary[label] = summary.get(label, 0) + 1

    return {
        "available": True,
        "uniprot": acc_used,
        "gene": data["gene"],
        "uniprot_length": len(up_seq),
        "variant_count": len(variants),
        "mapped_variant_count": mapped_count,
        "unmapped_variant_count": unmapped,
        "mapped_residue_count": len(residue_list),
        "residues": residue_list,
        "pathogenicity_summary": summary,
        "source": "EBI Proteins API (ClinVar clinical significance + gnomAD "
                  "allele frequencies)",
        "retrieved_utc": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%d %H:%M UTC"
        ),
        "locus": locus_info(data["gene"]),
    }


def locus_info(gene: str | None) -> dict | None:
    """Best-effort gene/locus card (chromosome, cytogenetic band) via NCBI Gene.

    Info-only and degrades gracefully — any failure returns None so the module
    never crashes on a locus lookup.
    """
    if not gene:
        return None
    try:
        import urllib.parse
        term = urllib.parse.quote(f"{gene}[sym] AND 9606[taxid]")
        search = fetch_json(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gene&term={term}&retmode=json"
        )
        ids = ((search.get("esearchresult") or {}).get("idlist")) or []
        if not ids:
            return {"gene": gene}
        gid = ids[0]
        summ = fetch_json(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gene&id={gid}&retmode=json"
        )
        doc = (summ.get("result") or {}).get(gid) or {}
        return {
            "gene": doc.get("name") or gene,
            "description": doc.get("description"),
            "chromosome": doc.get("chromosome"),
            "cytogenetic_band": doc.get("maplocation"),
            "ncbi_gene_id": gid,
        }
    except (FetchError, KeyError, ValueError, TypeError):
        return {"gene": gene}
