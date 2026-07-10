"""Tests for variant parsing + residue-level analysis (Phase 3). Offline."""

import unittest

from snaclex import variants
from snaclex.variants import VariantError


class TestParse(unittest.TestCase):
    def test_simple_one_letter(self):
        v = variants.parse("BRAF V600E")
        self.assertEqual(v["gene_or_acc"], "BRAF")
        self.assertEqual((v["wt"], v["position"], v["mut"]), ("V", 600, "E"))
        self.assertEqual(v["notation"], "V600E")

    def test_simple_three_letter(self):
        v = variants.parse("TP53 Arg175His")
        self.assertEqual((v["wt"], v["position"], v["mut"]), ("R", 175, "H"))

    def test_hgvs(self):
        v = variants.parse("P15056:p.Val600Glu")
        self.assertEqual(v["gene_or_acc"], "P15056")
        self.assertEqual(v["notation"], "V600E")

    def test_nonsense(self):
        v = variants.parse("BRCA1 Q1200*")
        self.assertEqual(v["mut"], "*")

    def test_synonymous_hgvs(self):
        v = variants.parse("P15056:p.Val600=")
        self.assertEqual(v["wt"], v["mut"])

    def test_bad_input(self):
        with self.assertRaises(VariantError):
            variants.parse("not a variant")

    def test_bad_amino_acid(self):
        with self.assertRaises(VariantError):
            variants.parse("BRAF Xyz600Glu")


class TestConsequence(unittest.TestCase):
    def test_missense(self):
        self.assertEqual(variants.consequence("V", "E"), "missense")

    def test_nonsense(self):
        self.assertIn("nonsense", variants.consequence("Q", "*"))

    def test_synonymous(self):
        self.assertEqual(variants.consequence("V", "V"), "synonymous")


def _record(seq, **over):
    rec = {
        "sequence": {"value": seq}, "review_status": "reviewed",
        "domains_motifs": [{"name": "Protein kinase", "type": "domain",
                            "start": 5, "end": 15}],
        "residue_annotations": [{"kind": "active_site", "start": 10, "end": 10,
                                 "description": "Proton acceptor"}],
        "ptms": [], "variants": [{"position": 10, "wt": "K", "mut": "R",
                                  "description": "in a disease; dbSNP"}],
        "disease_associations": [{"name": "Example syndrome",
                                  "description": "curated disease link"}],
        "structures": {"experimental": [{"pdb_id": "1ABC"}]},
    }
    rec.update(over)
    return rec


class TestAnalyze(unittest.TestCase):
    SEQ = "MDLSAKLICE" + "FGHIKLMNPQ"  # length 20; residue 10 = 'E'? check below

    def test_wt_match_and_domain(self):
        seq = "MDLSAKLICE"  # positions 1..10
        rec = _record(seq)
        # position 6 = 'K' (M1 D2 L3 S4 A5 K6 ...)
        v = {"wt": "K", "position": 6, "mut": "R", "notation": "K6R", "input": "X K6R"}
        out = variants.analyze(rec, v)
        self.assertTrue(out["sequence_mapping"]["wt_matches"])
        self.assertTrue(out["domain_disruption"])  # 6 within domain 5-15
        self.assertEqual(out["consequence"], "missense")

    def test_wt_mismatch_warns(self):
        seq = "MDLSAKLICE"
        rec = _record(seq)
        v = {"wt": "W", "position": 6, "mut": "R", "notation": "W6R", "input": "X W6R"}
        out = variants.analyze(rec, v)
        self.assertFalse(out["sequence_mapping"]["wt_matches"])
        self.assertTrue(any("mismatch" in w for w in out["warnings"]))

    def test_out_of_range(self):
        rec = _record("MDLSAKLICE")
        v = {"wt": "K", "position": 999, "mut": "R", "notation": "K999R", "input": "x"}
        out = variants.analyze(rec, v)
        self.assertFalse(out["sequence_mapping"]["wt_matches"])
        self.assertTrue(any("outside" in w for w in out["warnings"]))

    def test_residue_context_and_known_variant(self):
        rec = _record("MDLSAKLICE")
        v = {"wt": "I", "position": 10, "mut": "T", "notation": "I10T", "input": "x"}
        # seq[9] = 'E' actually; wt mismatch expected but we still get context at pos 10
        out = variants.analyze(rec, v)
        self.assertTrue(out["residue_context"])          # active_site at 10
        self.assertIsNotNone(out["known_variant"])       # known variant at 10

    def test_clinical_carries_evidence_and_disclaimer(self):
        rec = _record("MDLSAKLICE")
        v = {"wt": "E", "position": 10, "mut": "K", "notation": "E10K", "input": "x"}
        out = variants.analyze(rec, v)
        self.assertEqual(out["clinical"]["review_status"], "reviewed")
        self.assertTrue(out["clinical"]["interpretations"])
        self.assertIn("Not a clinical interpretation", out["clinical"]["disclaimer"])

    def test_structure_mapping_listed(self):
        rec = _record("MDLSAKLICE")
        v = {"wt": "K", "position": 6, "mut": "R", "notation": "K6R", "input": "x"}
        out = variants.analyze(rec, v)
        self.assertEqual(out["structure_mapping"]["experimental_structures"], ["1ABC"])


if __name__ == "__main__":
    unittest.main()
