"""Offline tests for homology-based structure selection (Phase 2)."""

import unittest
from unittest import mock

from snaclex import homology


class TestClassifyHits(unittest.TestCase):
    HITS = [
        {"pdb_id": "1IEP", "entity": "1", "identity": 0.99, "evalue": 0.0},
        {"pdb_id": "2GQG", "entity": "1", "identity": 0.62, "evalue": 1e-40},
        {"pdb_id": "3XYZ", "entity": "1", "identity": 0.25, "evalue": 0.1},
    ]

    def test_identity_preserved_and_gated(self):
        out = homology.classify_hits(self.HITS)
        self.assertEqual(len(out), 3)
        by_id = {h["pdb_id"]: h for h in out}
        self.assertTrue(by_id["1IEP"]["near_identical"])
        self.assertTrue(by_id["2GQG"]["docking_suitable"])       # 0.62 ≥ 0.40
        self.assertFalse(by_id["3XYZ"]["docking_suitable"])      # 0.25 < 0.40
        self.assertEqual(by_id["2GQG"]["sequence_identity"], 0.62)

    def test_direct_structures_excluded(self):
        out = homology.classify_hits(self.HITS, direct_pdb_ids=["1IEP"])
        self.assertNotIn("1IEP", {h["pdb_id"] for h in out})

    def test_deduplicates_on_pdb(self):
        dup = self.HITS + [{"pdb_id": "2gqg", "entity": "2", "identity": 0.62}]
        out = homology.classify_hits(dup)
        self.assertEqual(sum(1 for h in out if h["pdb_id"] == "2GQG"), 1)

    def test_max_report(self):
        many = [{"pdb_id": f"P{i:03d}", "identity": 0.5} for i in range(20)]
        self.assertEqual(len(homology.classify_hits(many, max_report=5)), 5)


class TestFindHomologousStructures(unittest.TestCase):
    def test_short_sequence_unavailable(self):
        out = homology.find_homologous_structures("ACDEFG")
        self.assertFalse(out["available"])

    def test_composes_search_and_classify(self):
        with mock.patch.object(homology.rcsb, "sequence_search",
                               return_value=[{"pdb_id": "2GQG", "identity": 0.62}]):
            out = homology.find_homologous_structures("A" * 40)
        self.assertTrue(out["available"])
        self.assertEqual(out["structures"][0]["pdb_id"], "2GQG")
        self.assertIn("identity_cutoff", out["method"])


if __name__ == "__main__":
    unittest.main()
