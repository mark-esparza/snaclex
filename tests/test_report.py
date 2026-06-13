"""Tests for the plain-language report summarizer."""

import unittest

from snaclex import report


def _profile(counts):
    return {
        "counts": counts,
        "component": {"label": "LIG B900", "kind": "ligand"},
        "contact_residues": [
            {"res_name": "ASP", "res_seq": 20},
            {"res_name": "PHE", "res_seq": 50},
        ],
        "contact_residue_count": 2,
    }


class TestReport(unittest.TestCase):
    def test_summary_mentions_pdb_and_component(self):
        out = report.summarize(
            _profile({"hydrogen_bond": 3}),
            {"pdb_id": "1ABC", "title": "Test protein"},
        )
        self.assertIn("1ABC", out["summary"])
        self.assertIn("LIG B900", out["summary"])
        self.assertIn("Research-only", out["disclaimer"])

    def test_hypotheses_generated_per_interaction(self):
        counts = {
            "metal_coordination": 1,
            "salt_bridge": 1,
            "hydrogen_bond": 3,
            "hydrophobic": 4,
            "aromatic": 1,
        }
        out = report.summarize(_profile(counts), {"pdb_id": "1ABC"})
        # Each interaction class contributes at least one hypothesis line.
        self.assertGreaterEqual(len(out["hypotheses"]), 5)

    def test_hypotheses_fallback_when_empty(self):
        empty = {
            "counts": {t: 0 for t in
                       ("metal_coordination", "salt_bridge", "hydrogen_bond",
                        "hydrophobic", "aromatic")},
            "component": {"label": "LIG B900", "kind": "ligand"},
            "contact_residues": [],
            "contact_residue_count": 0,
        }
        out = report.summarize(empty)
        self.assertEqual(len(out["hypotheses"]), 1)  # the "no strong interactions" line

    def test_handles_missing_metadata(self):
        out = report.summarize(_profile({"hydrogen_bond": 1}))
        self.assertIn("summary", out)


if __name__ == "__main__":
    unittest.main()
