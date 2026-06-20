"""Tests for the RCSB client: ID normalization + offline search parsing."""

import unittest
from unittest import mock

from snaclex import rcsb
from snaclex.http_util import FetchError, RateLimitError


class TestNormalizePdbId(unittest.TestCase):
    def test_valid_uppercased_and_trimmed(self):
        self.assertEqual(rcsb.normalize_pdb_id("  1abc "), "1ABC")

    def test_invalid_length_raises(self):
        with self.assertRaises(FetchError):
            rcsb.normalize_pdb_id("zz")

    def test_non_alnum_raises(self):
        with self.assertRaises(FetchError):
            rcsb.normalize_pdb_id("1a!c")


class TestFetchStructure(unittest.TestCase):
    def test_falls_back_to_mmcif_on_pdb_404(self):
        def fake(url, **kw):
            if url.endswith(".pdb"):
                raise FetchError("HTTP 404 for ...")
            return "data_6BCX\n_atom_site.group_PDB\n"
        with mock.patch.object(rcsb, "fetch_text", side_effect=fake):
            text = rcsb.fetch_structure("6BCX")
        self.assertIn("_atom_site.", text)

    def test_uses_pdb_when_available(self):
        with mock.patch.object(rcsb, "fetch_text", return_value="ATOM ...") as m:
            self.assertEqual(rcsb.fetch_structure("1HSG"), "ATOM ...")
            self.assertEqual(m.call_count, 1)  # no mmCIF fallback needed

    def test_rate_limit_does_not_trigger_fallback(self):
        with mock.patch.object(rcsb, "fetch_text",
                               side_effect=RateLimitError("throttled")) as m:
            with self.assertRaises(RateLimitError):
                rcsb.fetch_structure("1HSG")
            self.assertEqual(m.call_count, 1)  # surfaced, not retried as .cif


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


class TestFetchChemComponents(unittest.TestCase):
    def test_parses_name_formula_smiles_synonyms(self):
        fake = {"data": {"chem_comps": [
            {"chem_comp": {"id": "STI", "name": "Imatinib", "formula": "C29 H31 N7 O",
                           "formula_weight": 493.6, "type": "non-polymer",
                           "pdbx_synonyms": "Gleevec; Glivec"},
             "pdbx_chem_comp_descriptor": [
                 {"type": "InChI", "descriptor": "InChI=1S/..."},
                 {"type": "SMILES_CANONICAL", "descriptor": "Cc1ccc(cc1Nc1..."}]},
        ]}}
        with mock.patch.object(rcsb, "fetch_json", return_value=fake):
            out = rcsb.fetch_chem_components(["STI", "sti"])  # dedup/uppercase
        self.assertIn("STI", out)
        self.assertEqual(out["STI"]["name"], "Imatinib")
        self.assertEqual(out["STI"]["formula"], "C29 H31 N7 O")
        self.assertEqual(out["STI"]["smiles"], "Cc1ccc(cc1Nc1...")  # first SMILES
        self.assertEqual(out["STI"]["synonyms"], ["Gleevec", "Glivec"])

    def test_empty_for_no_codes(self):
        self.assertEqual(rcsb.fetch_chem_components([]), {})
        self.assertEqual(rcsb.fetch_chem_components(["", "  "]), {})

    def test_graceful_on_fetch_error(self):
        with mock.patch.object(rcsb, "fetch_json", side_effect=FetchError("down")):
            self.assertEqual(rcsb.fetch_chem_components(["ZN"]), {})

    def test_rest_fallback_when_graphql_empty(self):
        # GraphQL returns nothing; the per-code REST endpoint resolves the name.
        def fake(url, **kw):
            if "graphql" in url:
                return {"data": {"chem_comps": []}}
            return {"chem_comp": {"id": "ZN", "name": "ZINC ION", "formula": "Zn 2",
                                  "formula_weight": 65.4, "type": "non-polymer"},
                    "rcsb_chem_comp_synonyms": [{"name": "Zn2+"}],
                    "pdbx_chem_comp_descriptor": [{"type": "SMILES", "descriptor": "[Zn+2]"}]}
        with mock.patch.object(rcsb, "fetch_json", side_effect=fake):
            out = rcsb.fetch_chem_components(["ZN"])
        self.assertEqual(out["ZN"]["name"], "ZINC ION")
        self.assertEqual(out["ZN"]["smiles"], "[Zn+2]")
        self.assertEqual(out["ZN"]["synonyms"], ["Zn2+"])


if __name__ == "__main__":
    unittest.main()
