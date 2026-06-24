"""Tests for the Ensembl VEP genomic->protein consequence client (offline)."""

import os
import unittest
from unittest import mock

from snaclex import vep
from snaclex.http_util import FetchError


class TestAvailability(unittest.TestCase):
    def test_off_by_default(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertFalse(vep.available())

    def test_on_when_enabled(self):
        with mock.patch.dict(os.environ, {"SNACLEX_ENABLE_VEP": "1"}, clear=True):
            self.assertTrue(vep.available())
        with mock.patch.dict(os.environ, {"SNACLEX_ENABLE_VEP": "off"}, clear=True):
            self.assertFalse(vep.available())


# A trimmed VEP region response for chr17:7676154 G>A (TP53-like).
_PAYLOAD = [{
    "transcript_consequences": [
        {"gene_symbol": "OTHER", "consequence_terms": ["intron_variant"]},
        {"gene_symbol": "TP53", "canonical": 1, "swissprot": ["P04637.4"],
         "protein_start": 273, "amino_acids": "R/H",
         "consequence_terms": ["missense_variant"]},
    ]
}]


class TestPickConsequence(unittest.TestCase):
    def test_picks_canonical_protein_change(self):
        out = vep._pick_consequence(_PAYLOAD)
        self.assertEqual(out["gene"], "TP53")
        self.assertEqual(out["uniprot"], "P04637")   # version stripped
        self.assertEqual(out["protein_position"], 273)
        self.assertEqual((out["wt_aa"], out["mut_aa"]), ("R", "H"))
        self.assertEqual(out["consequence"], "missense_variant")

    def test_none_when_no_protein_change(self):
        payload = [{"transcript_consequences": [
            {"gene_symbol": "X", "consequence_terms": ["intron_variant"]}]}]
        self.assertIsNone(vep._pick_consequence(payload))

    def test_none_for_synonymous(self):
        payload = [{"transcript_consequences": [
            {"swissprot": "P1", "protein_start": 5, "amino_acids": "L/L"}]}]
        self.assertIsNone(vep._pick_consequence(payload))


class TestAnnotateOne(unittest.TestCase):
    def test_success_with_injected_fetcher(self):
        out = vep.annotate_one("17", 7676154, "G", "A", fetcher=lambda url: _PAYLOAD)
        self.assertEqual(out["uniprot"], "P04637")

    def test_graceful_on_fetch_error(self):
        def boom(url):
            raise FetchError("blocked")
        self.assertIsNone(vep.annotate_one("1", 100, "A", "T", fetcher=boom))

    def test_region_url_shape(self):
        url = vep._region_url("chr17", 7676154, "G", "A")
        self.assertIn("/vep/human/region/17:7676154-7676154/A", url)


if __name__ == "__main__":
    unittest.main()
