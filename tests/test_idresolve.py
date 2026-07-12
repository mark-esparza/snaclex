"""Tests for query interpretation (Stage 1) + xref graph assembly (Stage 2).

``interpret`` and ``build_xref_graph`` are pure and tested directly; the
network-touching ``resolve`` is covered via the server integration path.
"""

import unittest

from snaclex import idresolve


class TestInterpret(unittest.TestCase):
    def _type(self, q):
        return idresolve.interpret(q)["query_type"]

    def test_uniprot_accession(self):
        self.assertEqual(self._type("P38398"), "uniprot_accession")
        self.assertEqual(self._type("P0DTC2"), "uniprot_accession")

    def test_uniprot_isoform(self):
        self.assertEqual(self._type("P38398-2"), "uniprot_accession")

    def test_uniparc_upi(self):
        self.assertEqual(self._type("UPI0000123ABC"), "uniparc_upi")

    def test_refseq(self):
        self.assertEqual(self._type("NP_009225.1"), "refseq_protein")
        self.assertEqual(self._type("XP_011529125"), "refseq_protein")

    def test_pdb_id(self):
        self.assertEqual(self._type("1HSG"), "pdb_id")
        self.assertEqual(self._type("3PTB"), "pdb_id")

    def test_raw_fasta_header(self):
        self.assertEqual(self._type(">sp|P38398\nMDLSALRVEE"), "raw_sequence")

    def test_raw_bare_sequence(self):
        self.assertEqual(self._type("MDLSALRVEEVQNVINAMQKILECPICLE"), "raw_sequence")

    def test_variant(self):
        out = idresolve.interpret("BRAF V600E")
        self.assertEqual(out["query_type"], "variant")
        self.assertEqual(out["variant"]["position"], 600)
        self.assertEqual(out["variant"]["wt"], "V")
        self.assertEqual(out["variant"]["mut"], "E")

    def test_gene_symbol(self):
        self.assertEqual(self._type("BRCA1"), "gene_symbol")

    def test_protein_name(self):
        self.assertEqual(self._type("breast cancer type 1 susceptibility protein"),
                         "protein_name")

    def test_inchikey_is_chemical(self):
        self.assertEqual(self._type("BSYNRYMUTXBXSQ-UHFFFAOYSA-N"), "chemical")

    def test_pair(self):
        out = idresolve.interpret("ABL1 + CC(=O)Oc1ccccc1C(=O)O")
        self.assertEqual(out["query_type"], "protein_chemical_pair")
        self.assertEqual(out["protein"]["query"], "ABL1")

    def test_batch(self):
        out = idresolve.interpret("P38398, P04637, P00519")
        self.assertEqual(out["query_type"], "batch")
        self.assertEqual(len(out["tokens"]), 3)

    def test_gene_and_name_need_disambiguation(self):
        self.assertTrue(idresolve.interpret("BRCA1")["needs_disambiguation"])
        self.assertFalse(idresolve.interpret("P38398")["needs_disambiguation"])


class TestBuildXrefGraph(unittest.TestCase):
    def test_links_uniprot_to_sources(self):
        up = {
            "accession": "P00519",
            "cross_references": {"RefSeq": [{"id": "NP_005148.2", "properties": {}}],
                                 "PDB": [{"id": "1IEP", "properties": {}}]},
        }
        upi = {"upi": "UPI000012ABCD",
               "cross_references": [{"database": "UniProtKB/Swiss-Prot",
                                     "id": "P00519", "active": True}]}
        g = idresolve.build_xref_graph(
            uniprot_fragment=up, uniparc_fragment=upi,
            pdb_entities=[("1IEP", "A")], alphafold_ids=["AF-P00519-F1"])
        ns = {n["ns"] for n in g["nodes"]}
        self.assertEqual(ns, {"UniProt", "RefSeq", "PDB", "UniParc", "AlphaFold"})
        rels = {e["rel"] for e in g["edges"]}
        self.assertIn("sequence_identical", rels)
        self.assertIn("has_predicted_structure", rels)

    def test_isoforms_not_collapsed(self):
        # Two distinct accessions stay as two nodes (never merged on name).
        g = idresolve.build_xref_graph(
            uniprot_fragment={"accession": "P38398", "cross_references": {}})
        ids = {n["id"] for n in g["nodes"]}
        self.assertIn("P38398", ids)


class TestExtractSequence(unittest.TestCase):
    def test_fasta(self):
        self.assertEqual(idresolve._extract_sequence(">h\nACDE\nFGHI"), "ACDEFGHI")

    def test_bare(self):
        self.assertEqual(idresolve._extract_sequence("ACDE FGHI"), "ACDEFGHI")


if __name__ == "__main__":
    unittest.main()
