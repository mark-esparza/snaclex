"""Tests for the CELL×GENE expression-context client (offline)."""

import os
import unittest
from unittest import mock

from snaclex import cellxgene
from snaclex.http_util import FetchError


class TestLinkAndAvailability(unittest.TestCase):
    def test_link_always_built(self):
        link = cellxgene.discover_link("HLA-A")
        self.assertIn("gene-expression?genes=HLA-A", link)

    def test_off_by_default(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertFalse(cellxgene.available())
            out = cellxgene.gene_expression("HLA-A")
            self.assertFalse(out["available"])
            self.assertIn("gene-expression", out["link"])  # link still offered

    def test_on_when_enabled(self):
        with mock.patch.dict(os.environ, {"SNACLEX_ENABLE_CELLXGENE": "1"}, clear=True):
            self.assertTrue(cellxgene.available())


class TestGeneExpression(unittest.TestCase):
    def _fetcher(self, xref, wmg):
        def f(url, body=None):
            if "xrefs/symbol" in url:
                return xref
            return wmg
        return f

    def test_success(self):
        xref = [{"id": "ENSG00000206503", "type": "gene"}]
        wmg = {
            "expression_summary": {"ENSG00000206503": {
                "UBERON:0000178": {  # blood
                    "CL:0000236": {"me": 3.1, "pc": 0.8},  # B cell
                    "CL:0000084": {"me": 2.2, "pc": 0.5},  # T cell
                }}},
            "term_id_labels": {"cell_types": {
                "CL:0000236": {"name": "B cell"}, "CL:0000084": {"name": "T cell"}}},
        }
        out = cellxgene.gene_expression("HLA-A", fetcher=self._fetcher(xref, wmg))
        self.assertTrue(out["available"])
        self.assertEqual(out["ensembl_id"], "ENSG00000206503")
        self.assertEqual(out["cell_types"][0]["cell_type"], "B cell")  # highest expr
        self.assertEqual(out["cell_types"][0]["mean_expr"], 3.1)

    def test_unresolved_gene(self):
        out = cellxgene.gene_expression("NOTAGENE", fetcher=self._fetcher([], {}))
        self.assertFalse(out["available"])
        self.assertIn("Ensembl", out["reason"])

    def test_graceful_on_wmg_error(self):
        def f(url, body=None):
            if "xrefs/symbol" in url:
                return [{"id": "ENSG1"}]
            raise FetchError("blocked")
        out = cellxgene.gene_expression("HLA-A", fetcher=f)
        self.assertFalse(out["available"])
        self.assertIn("unavailable", out["reason"])
        self.assertIn("gene-expression", out["link"])

    def test_empty_summary(self):
        out = cellxgene.gene_expression("HLA-A", fetcher=self._fetcher([{"id": "ENSG1"}], {}))
        self.assertFalse(out["available"])


if __name__ == "__main__":
    unittest.main()
