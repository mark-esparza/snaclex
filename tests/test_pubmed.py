"""Offline tests for the PubMed literature adapter + UniProt reference extract."""

import unittest

from snaclex import pubmed, uniprot


class TestParseElink(unittest.TestCase):
    def test_extracts_pmids(self):
        data = {"linksets": [{"linksetdbs": [
            {"linkname": "protein_pubmed", "links": ["111", "222", "222"]}]}]}
        self.assertEqual(pubmed.parse_elink(data), ["111", "222"])

    def test_empty(self):
        self.assertEqual(pubmed.parse_elink({}), [])


class TestParseEsummary(unittest.TestCase):
    DATA = {"result": {"uids": ["8552191"], "8552191": {
        "title": "A landmark paper.",
        "source": "Nature",
        "pubdate": "1994 Jan 1",
        "authors": [{"name": "Smith J"}, {"name": "Doe A"}, {"name": "Roe B"},
                    {"name": "Poe C"}],
        "articleids": [{"idtype": "pubmed", "value": "8552191"},
                       {"idtype": "doi", "value": "10.1000/x"}],
    }}}

    def test_fields(self):
        out = pubmed.parse_esummary(self.DATA)
        self.assertEqual(len(out), 1)
        r = out[0]
        self.assertEqual(r["pmid"], "8552191")
        self.assertEqual(r["journal"], "Nature")
        self.assertEqual(r["year"], "1994")
        self.assertEqual(r["doi"], "10.1000/x")
        self.assertTrue(r["authors"].endswith("et al."))
        self.assertIn("8552191", r["url"])


class TestUniprotReferences(unittest.TestCase):
    ENTRY = {
        "primaryAccession": "P00519",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "sequence": {"value": "MLE", "length": 3},
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "references": [
            {"citation": {"title": "Cloning of ABL1.", "journal": "Cell",
                          "publicationDate": "1986",
                          "authors": ["Shtivelman E", "Lifshitz B"],
                          "citationCrossReferences": [
                              {"database": "PubMed", "id": "3018722"},
                              {"database": "DOI", "id": "10.1016/x"}]}},
            {"citation": {"title": "No-PMID ref.", "journal": "J",
                          "authors": ["X Y"], "citationCrossReferences": []}},
        ],
    }

    def test_literature_extracted_into_record(self):
        frag = uniprot.parse_entry(self.ENTRY)["data"]
        lit = frag["literature"]
        self.assertEqual(len(lit), 2)
        self.assertEqual(lit[0]["pmid"], "3018722")
        self.assertEqual(lit[0]["doi"], "10.1016/x")
        self.assertEqual(lit[0]["journal"], "Cell")
        self.assertEqual(lit[0]["source"], "UniProtKB reference")

    def test_dedup_by_pmid(self):
        entry = dict(self.ENTRY)
        entry["references"] = self.ENTRY["references"] + [self.ENTRY["references"][0]]
        lit = uniprot.parse_entry(entry)["data"]["literature"]
        pmids = [x["pmid"] for x in lit if x["pmid"]]
        self.assertEqual(len(pmids), len(set(pmids)))


if __name__ == "__main__":
    unittest.main()
