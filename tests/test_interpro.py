"""Offline tests for the InterPro enrichment adapter (Phase 2)."""

import os
import unittest
from unittest import mock

from snaclex import interpro


class TestEnabled(unittest.TestCase):
    def test_disabled_by_default(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertFalse(interpro.enabled())

    def test_enabled_flag(self):
        for val in ("1", "true", "on", "YES"):
            with mock.patch.dict(os.environ, {"SNACLEX_ENABLE_INTERPRO": val}):
                self.assertTrue(interpro.enabled())


class TestFetchGraceful(unittest.TestCase):
    def test_disabled_returns_unavailable(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            out = interpro.fetch_domains("P00519")
        self.assertFalse(out["available"])
        self.assertEqual(out["entries"], [])
        self.assertIn("disabled", out["reason"])


class TestParseEntries(unittest.TestCase):
    DATA = {
        "results": [
            {"metadata": {"accession": "IPR000719", "name": "Protein kinase domain",
                          "type": "domain", "source_database": "interpro"},
             "proteins": [{"entry_protein_locations": [
                 {"fragments": [{"start": 242, "end": 493}]}]}]},
            {"metadata": {"accession": "PF07714", "name": "PK_Tyr_Ser-Thr",
                          "type": "domain", "source_database": "pfam"},
             "proteins": [{"entry_protein_locations": [
                 {"fragments": [{"start": 248, "end": 500}]}]}]},
        ]
    }

    def test_parses_entries_and_locations(self):
        entries = interpro.parse_entries(self.DATA)
        self.assertEqual(len(entries), 2)
        e0 = entries[0]
        self.assertEqual(e0["accession"], "IPR000719")
        self.assertEqual(e0["type"], "domain")
        self.assertEqual(e0["locations"], [{"start": 242, "end": 493}])
        self.assertEqual(e0["source"], "InterPro")

    def test_empty_data(self):
        self.assertEqual(interpro.parse_entries({}), [])


if __name__ == "__main__":
    unittest.main()
