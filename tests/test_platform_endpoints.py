"""Integration tests for the sequence-first platform endpoints.

Boots the real in-process server and mocks the network-touching resolver /
structure-availability / PubChem calls, so the endpoint composition (routing →
build_record → evidence separation) is exercised end-to-end without any live
upstream call — matching tests/test_server_integration.py conventions.
"""

import http.client
import json
import threading
import unittest
from http.server import ThreadingHTTPServer
from unittest import mock

import server
from snaclex import uniprot

_ENTRY = {
    "primaryAccession": "P00519", "uniProtkbId": "ABL1_HUMAN",
    "entryType": "UniProtKB reviewed (Swiss-Prot)", "entryAudit": {"entryVersion": 1},
    "sequence": {"value": "MLEICLKLVGCKSKKGLSSS", "length": 20, "crc64": "AABB"},
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "genes": [{"geneName": {"value": "ABL1"}}],
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Tyrosine-protein kinase ABL1"}}},
    "comments": [], "features": [],
    "uniProtKBCrossReferences": [{"database": "GeneID", "id": "25", "properties": []}],
}


def _resolution_with_uniprot():
    return {
        "query": "P00519", "query_type": "uniprot_accession",
        "chosen": {"accession": "P00519", "taxon_id": 9606, "review_status": "reviewed"},
        "xref_graph": {"nodes": [], "edges": []}, "identity_key": {"crc64": "AABB"},
        "warnings": [], "_uniprot_fragment": uniprot.parse_entry(_ENTRY)["data"],
        "_uniparc_fragment": None,
    }


class TestPlatformEndpoints(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.port = cls.httpd.server_address[1]
        cls.thread = threading.Thread(target=cls.httpd.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()

    def setUp(self):
        server._IP_LIMITER._buckets.clear()
        server._HEAVY_LIMITER._buckets.clear()

    def _get(self, path):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("GET", path)
        resp = conn.getresponse()
        body = resp.read()
        conn.close()
        return resp, json.loads(body)

    def test_protein_experimental_structure(self):
        avail = {"experimental": [{"pdb_id": "1IEP", "kind": "experimental",
                                   "docking_suitable": True}],
                 "homologous": [], "predicted": [], "tier": "experimental",
                 "best_for_docking": "1IEP", "warnings": []}
        with mock.patch.object(server.idresolve, "resolve",
                               return_value=_resolution_with_uniprot()), \
             mock.patch.object(server.proteinrecord, "structure_availability",
                               return_value=avail):
            resp, rec = self._get("/api/protein?q=P00519")
        self.assertEqual(resp.status, 200)
        self.assertEqual(rec["canonical_id"], "SNX:PRT:9606:P00519:1")
        self.assertEqual(rec["structures"]["tier"], "experimental")
        self.assertIsNotNone(rec["sequence_analysis"])  # always present

    def test_protein_alphafold_fallback(self):
        avail = {"experimental": [], "homologous": [],
                 "predicted": [{"model_id": "AF-P00519-F1", "docking_suitable": True}],
                 "tier": "predicted", "best_for_docking": "AF-P00519-F1",
                 "warnings": ["structure is a predicted model"]}
        with mock.patch.object(server.idresolve, "resolve",
                               return_value=_resolution_with_uniprot()), \
             mock.patch.object(server.proteinrecord, "structure_availability",
                               return_value=avail):
            resp, rec = self._get("/api/protein?q=P00519")
        self.assertEqual(rec["structures"]["tier"], "predicted")
        self.assertTrue(rec["structures"]["predicted"])

    def test_protein_sequence_only(self):
        res = {"query": ">x\nMLEICLKLVGCKSKKGLSSS", "query_type": "raw_sequence",
               "sequence": "MLEICLKLVGCKSKKGLSSS", "identity_key": {"crc64": "AABB"},
               "warnings": ["novel/unmatched sequence — sequence-only"],
               "xref_graph": {"nodes": [], "edges": []}}
        with mock.patch.object(server.idresolve, "resolve", return_value=res):
            resp, rec = self._get("/api/protein?q=RAWSEQ")
        self.assertEqual(resp.status, 200)
        self.assertIsNotNone(rec["sequence_analysis"])
        self.assertEqual(rec["structures"]["tier"], "sequence_only")

    def test_sequence_analysis_endpoint(self):
        with mock.patch.object(server.idresolve, "resolve",
                               return_value=_resolution_with_uniprot()):
            resp, out = self._get("/api/sequence_analysis?acc=P00519&ph=7.4")
        self.assertEqual(resp.status, 200)
        self.assertEqual(out["evidence_category"], "calculated")
        self.assertEqual(out["sequence_analysis"]["charge_at_ph"]["ph"], 7.4)

    def test_evidence_separates_levels(self):
        assay_rows = [
            {"aid": "372", "outcome": "Active", "target_geneid": "25",
             "activity_name": "IC50", "activity_value_uM": 0.038,
             "assay_type": "confirmatory", "pmid": "12345"},
            {"aid": "999", "outcome": "Inactive", "target_geneid": "25",
             "activity_name": "IC50", "activity_value_uM": None,
             "assay_type": "screening", "pmid": None},
        ]
        avail = {"tier": "experimental", "best_for_docking": "1IEP",
                 "experimental": [], "predicted": [], "warnings": []}
        with mock.patch.object(server.idresolve, "resolve",
                               return_value=_resolution_with_uniprot()), \
             mock.patch.object(server.pubchem, "lookup_compound",
                               return_value={"cid": 5291}), \
             mock.patch.object(server.pubchem, "fetch_identity",
                               return_value={"cid": 5291, "inchikey": "KTUFNOKKBVMGRW-UHFFFAOYSA-N"}), \
             mock.patch.object(server.pubchem, "fetch_parent_cid", return_value=None), \
             mock.patch.object(server.pubchem, "bioassay_summary", return_value=assay_rows), \
             mock.patch.object(server.proteinrecord, "structure_availability",
                               return_value=avail):
            resp, out = self._get("/api/evidence?protein=P00519&chemical=imatinib")
        self.assertEqual(resp.status, 200)
        # Two Level-B rows (including the inactive one) attributed by GeneID.
        self.assertTrue(all(e["level"] == "B" for e in out["evidence"]))
        self.assertEqual(out["summary"]["by_level"]["B"], 2)
        self.assertNotIn("interacts", out["summary"])  # never a boolean
        # Docking is a SEPARATE hypothesis, not merged into the evidence list.
        self.assertTrue(out["docking"]["available"])
        self.assertIn("fit score", out["docking"]["note"])

    def test_evidence_requires_both_params(self):
        resp, out = self._get("/api/evidence?protein=P00519")
        self.assertEqual(resp.status, 400)

    def test_domains_separates_sources(self):
        res = _resolution_with_uniprot()
        res["_uniprot_fragment"]["domains_motifs"] = [
            {"type": "domain", "name": "Protein kinase", "start": 242, "end": 493}]
        interpro_out = {"available": True, "entries": [
            {"accession": "IPR000719", "name": "Protein kinase domain",
             "type": "domain", "source": "InterPro", "locations": [{"start": 242, "end": 493}]}]}
        with mock.patch.object(server.idresolve, "resolve", return_value=res), \
             mock.patch.object(server.interpro, "fetch_domains", return_value=interpro_out):
            resp, out = self._get("/api/domains?acc=P00519")
        self.assertEqual(resp.status, 200)
        self.assertEqual(out["curated"]["source"], "UniProtKB")
        self.assertEqual(out["curated"]["domains"][0]["source"], "UniProtKB")
        self.assertTrue(out["interpro"]["available"])
        self.assertEqual(out["interpro"]["entries"][0]["source"], "InterPro")

    def test_model_confidence_gates_docking(self):
        models = [{"data": {"model_id": "AF-P00519-F1",
                            "pdb_url": "https://x/AF-P00519-F1-model_v4.pdb",
                            "version": 4, "coverage_fraction": 1.0,
                            "fragmented": False, "pae_available": True},
                   "provenance": {"source": "AlphaFold DB"}}]

        def _ca(serial, b):
            line = list(" " * 80)

            def put(col, s):
                for i, ch in enumerate(s):
                    line[col - 1 + i] = ch
            put(1, "ATOM"); put(7, f"{serial:>5}"); put(13, "CA"); put(18, "ALA")
            put(22, "A"); put(23, f"{serial:>4}")
            put(31, f"{0.0:>8.3f}"); put(39, f"{0.0:>8.3f}"); put(47, f"{0.0:>8.3f}")
            put(55, f"{1.0:>6.2f}"); put(61, f"{b:>6.2f}"); put(77, " C")
            return "".join(line)

        pdb = "\n".join(_ca(i, 40.0) for i in range(1, 6))  # all low confidence
        with mock.patch.object(server.alphafold, "fetch_models", return_value=models), \
             mock.patch.object(server.alphafold, "fetch_model_pdb", return_value=pdb):
            resp, out = self._get("/api/model_confidence?acc=P00519")
        self.assertEqual(resp.status, 200)
        self.assertTrue(out["available"])
        self.assertFalse(out["confidence"]["docking_suitable"])
        self.assertTrue(out["confidence"]["low_confidence_regions"])

    def test_model_confidence_absent(self):
        with mock.patch.object(server.alphafold, "fetch_models", return_value=[]):
            resp, out = self._get("/api/model_confidence?acc=NOPE99")
        self.assertEqual(resp.status, 200)
        self.assertFalse(out["available"])


if __name__ == "__main__":
    unittest.main()
