"""Tests for the RCSB client: ID normalization + offline search parsing."""

import unittest
from unittest import mock

from snaclex import rcsb
from snaclex.http_util import FetchError


class TestNormalizePdbId(unittest.TestCase):
    def test_valid_uppercased_and_trimmed(self):
        self.assertEqual(rcsb.normalize_pdb_id("  1abc "), "1ABC")

    def test_invalid_length_raises(self):
        with self.assertRaises(FetchError):
            rcsb.normalize_pdb_id("zz")

    def test_non_alnum_raises(self):
        with self.assertRaises(FetchError):
            rcsb.normalize_pdb_id("1a!c")


class TestSearchByName(unittest.TestCase):
    def test_parses_result_set_even_if_enrichment_fails(self):
        fake = {"result_set": [
            {"identifier": "1ABC", "score": 1.0},
            {"identifier": "2XYZ", "score": 0.5},
        ]}
        with mock.patch.object(rcsb, "fetch_json", return_value=fake), \
             mock.patch.object(rcsb, "fetch_entry_summaries",
                               side_effect=FetchError("enrich down")):
            results = rcsb.search_by_name("kinase", limit=2)
        self.assertEqual([r["pdb_id"] for r in results], ["1ABC", "2XYZ"])
        self.assertEqual(results[0]["score"], 1.0)

    def test_enrichment_merged_when_available(self):
        fake = {"result_set": [{"identifier": "1ABC", "score": 1.0}]}
        summaries = {"1ABC": {"title": "A protein", "organism": "Homo sapiens"}}
        with mock.patch.object(rcsb, "fetch_json", return_value=fake), \
             mock.patch.object(rcsb, "fetch_entry_summaries", return_value=summaries):
            results = rcsb.search_by_name("kinase", limit=1)
        self.assertEqual(results[0]["title"], "A protein")
        self.assertEqual(results[0]["organism"], "Homo sapiens")


if __name__ == "__main__":
    unittest.main()
