"""Tests for the TOPMed/BRAVO client (network mocked at the fetch seam)."""

import unittest
from unittest import mock

from snaclex import bravo
from snaclex.http_util import FetchError


class TestCoerceFreq(unittest.TestCase):
    def test_direct_fields(self):
        out = bravo._coerce_freq({"allele_freq": 0.012, "allele_count": 30, "allele_num": 2500})
        self.assertAlmostEqual(out["allele_freq"], 0.012)

    def test_computed_from_counts(self):
        out = bravo._coerce_freq({"ac": 5, "an": 1000})
        self.assertAlmostEqual(out["allele_freq"], 0.005)

    def test_data_envelope(self):
        out = bravo._coerce_freq({"data": [{"allele_freq": 0.3}]})
        self.assertAlmostEqual(out["allele_freq"], 0.3)

    def test_missing_returns_none(self):
        self.assertIsNone(bravo._coerce_freq({"chrom": "1"}))


class TestClassifyFrequency(unittest.TestCase):
    def test_buckets(self):
        self.assertEqual(bravo.classify_frequency(None), "unknown")
        self.assertEqual(bravo.classify_frequency(0), "absent")
        self.assertIn("ultra-rare", bravo.classify_frequency(0.00001))
        self.assertIn("rare", bravo.classify_frequency(0.005))
        self.assertIn("common", bravo.classify_frequency(0.2))


class TestVariantFrequency(unittest.TestCase):
    def test_success(self):
        with mock.patch.object(bravo, "fetch_json", return_value={"allele_freq": 0.04}):
            out = bravo.variant_frequency("chr17", 7676154, "G", "A")
        self.assertTrue(out["available"])
        self.assertEqual(out["variant_id"], "17-7676154-G-A")
        self.assertAlmostEqual(out["allele_freq"], 0.04)

    def test_fetch_failure_is_graceful(self):
        with mock.patch.object(bravo, "fetch_json", side_effect=FetchError("403")):
            out = bravo.variant_frequency("1", 100, "A", "T")
        self.assertFalse(out["available"])
        self.assertIn("unavailable", out["reason"])

    def test_not_found_is_graceful(self):
        with mock.patch.object(bravo, "fetch_json", return_value={"chrom": "1"}):
            out = bravo.variant_frequency("1", 100, "A", "T")
        self.assertFalse(out["available"])


if __name__ == "__main__":
    unittest.main()
