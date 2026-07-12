"""Tests for protein comparison (family / ortholog views). Pure, offline."""

import unittest

from snaclex import compare


def _rec(acc, gene, taxon, organism, length, domains, tier="experimental"):
    return {
        "canonical_id": f"SNX:PRT:{taxon}:{acc}:1",
        "accessions": {"uniprot_primary": acc},
        "preferred_name": f"{gene} protein",
        "gene": {"symbol": gene},
        "organism": {"scientific_name": organism, "taxon_id": taxon},
        "review_status": "reviewed",
        "sequence": {"length": length},
        "sequence_analysis": {"length": length, "molecular_weight_Da": length * 110.0,
                              "isoelectric_point": 6.5, "gravy": -0.2},
        "domains_motifs": [{"name": d} for d in domains],
        "variants": [], "structures": {"tier": tier},
    }


HUMAN = _rec("P00519", "ABL1", 9606, "Homo sapiens", 1130, ["Protein kinase", "SH2"])
MOUSE = _rec("P00520", "ABL1", 10090, "Mus musculus", 1123, ["Protein kinase", "SH2"])
OTHER = _rec("P12931", "SRC", 9606, "Homo sapiens", 536, ["Protein kinase"])


class TestRecordSummary(unittest.TestCase):
    def test_summary_fields(self):
        s = compare.record_summary(HUMAN)
        self.assertEqual(s["accession"], "P00519")
        self.assertEqual(s["gene"], "ABL1")
        self.assertEqual(s["n_domains"], 2)
        self.assertEqual(s["structure_tier"], "experimental")


class TestSharedDomains(unittest.TestCase):
    def test_in_all_vs_some(self):
        out = compare.shared_domains([HUMAN, MOUSE, OTHER])
        self.assertIn("Protein kinase", out["in_all"])
        self.assertIn("SH2", out["in_some"])   # only ABL1 records have SH2
        self.assertNotIn("SH2", out["in_all"])


class TestOrthology(unittest.TestCase):
    def test_ortholog_candidates_cross_taxon(self):
        out = compare.orthology_groups([HUMAN, MOUSE, OTHER])
        self.assertIn("ABL1", out["ortholog_candidates"])       # human + mouse
        self.assertNotIn("SRC", out["ortholog_candidates"])     # single organism
        self.assertIn("candidate", out["note"].lower())

    def test_never_merges_on_symbol(self):
        # Grouping keeps distinct accessions/taxa, does not collapse them.
        out = compare.orthology_groups([HUMAN, MOUSE])
        members = out["ortholog_candidates"]["ABL1"]
        self.assertEqual({m["taxon_id"] for m in members}, {9606, 10090})


class TestBuildComparison(unittest.TestCase):
    def test_shape(self):
        out = compare.build_comparison([HUMAN, MOUSE, OTHER])
        self.assertEqual(out["n"], 3)
        self.assertEqual(len(out["rows"]), 3)
        self.assertEqual(out["length_range"], {"min": 536, "max": 1130})
        self.assertIn("columns", out)

    def test_empty(self):
        out = compare.build_comparison([])
        self.assertEqual(out["n"], 0)
        self.assertIsNone(out["length_range"])


if __name__ == "__main__":
    unittest.main()
