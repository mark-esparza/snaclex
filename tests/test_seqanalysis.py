"""Tests for the sequence-only analysis (Stage 3) — pure compute, offline."""

import unittest

from snaclex import seqanalysis


class TestCleanSequence(unittest.TestCase):
    def test_strips_whitespace_and_numbers(self):
        self.assertEqual(seqanalysis.clean_sequence(" ACD 12\nEFG "), "ACDEFG")

    def test_uppercases(self):
        self.assertEqual(seqanalysis.clean_sequence("acdefg"), "ACDEFG")

    def test_rejects_empty(self):
        with self.assertRaises(seqanalysis.SequenceError):
            seqanalysis.clean_sequence("   ")

    def test_rejects_non_protein(self):
        with self.assertRaises(seqanalysis.SequenceError):
            seqanalysis.clean_sequence("123456789")


class TestMolecularWeight(unittest.TestCase):
    def test_single_glycine(self):
        # G residue 57.0519 + water 18.01524 ≈ 75.07
        self.assertAlmostEqual(seqanalysis.molecular_weight("G"), 75.07, places=1)

    def test_scales_with_length(self):
        self.assertGreater(
            seqanalysis.molecular_weight("AAAA"), seqanalysis.molecular_weight("A"))


class TestIsoelectricPoint(unittest.TestCase):
    def test_neutral_backbone_only(self):
        # Only termini ionize → pI midway between C-term (3.6) and N-term (8.6).
        self.assertAlmostEqual(seqanalysis.isoelectric_point("AAAA"), 6.1, places=1)

    def test_basic_sequence_high_pi(self):
        self.assertGreater(seqanalysis.isoelectric_point("KKKKKK"), 9.5)

    def test_acidic_sequence_low_pi(self):
        self.assertLess(seqanalysis.isoelectric_point("DDDDDD"), 4.5)


class TestNetCharge(unittest.TestCase):
    def test_charge_zero_at_pi(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        pi = seqanalysis.isoelectric_point(seq)
        self.assertAlmostEqual(seqanalysis.net_charge(seq, pi), 0.0, places=1)

    def test_positive_below_pi(self):
        seq = "KKKKAAAA"
        pi = seqanalysis.isoelectric_point(seq)
        self.assertGreater(seqanalysis.net_charge(seq, pi - 2), 0)


class TestGravyAndProfile(unittest.TestCase):
    def test_gravy_hydrophobic_positive(self):
        self.assertGreater(seqanalysis.gravy("IIIIVVVV"), 0)

    def test_gravy_hydrophilic_negative(self):
        self.assertLess(seqanalysis.gravy("DDDDKKKK"), 0)

    def test_profile_length(self):
        prof = seqanalysis.hydropathy_profile("A" * 30, window=9)
        self.assertEqual(len(prof), 30 - 9 + 1)

    def test_profile_short_sequence_empty(self):
        self.assertEqual(seqanalysis.hydropathy_profile("AAA", window=9), [])


class TestLowComplexity(unittest.TestCase):
    def test_homopolymer_flagged(self):
        regions = seqanalysis.low_complexity_regions("Q" * 30, window=20)
        self.assertTrue(regions)
        self.assertEqual(regions[0]["start"], 1)

    def test_diverse_not_flagged(self):
        seq = ("ACDEFGHIKLMNPQRSTVWY" * 2)
        self.assertEqual(seqanalysis.low_complexity_regions(seq, window=20), [])


class TestAnalyze(unittest.TestCase):
    def test_full_shape(self):
        out = seqanalysis.analyze("ACDEFGHIKLMNPQRSTVWY", ph=7.4)
        self.assertEqual(out["length"], 20)
        self.assertIn("molecular_weight_Da", out)
        self.assertIn("isoelectric_point", out)
        self.assertEqual(out["charge_at_ph"]["ph"], 7.4)
        self.assertIn("method", out)
        self.assertIn("composition", out)

    def test_downsampling_caps_profile(self):
        out = seqanalysis.analyze("A" * 5000, profile_points=100)
        self.assertLessEqual(len(out["hydropathy_profile"]["downsampled"]), 101)


if __name__ == "__main__":
    unittest.main()
