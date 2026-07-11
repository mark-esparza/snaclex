"""Integration tests for the sequence-first platform endpoints.

Boots the real in-process server and mocks the network-touching resolver /
structure-availability / PubChem calls, so the endpoint composition (routing →
build_record → evidence separation) is exercised end-to-end without any live
upstream call — matching tests/test_server_integration.py conventions.
"""

import http.client
import json
import threading
import time
import unittest
from http.server import ThreadingHTTPServer
from unittest import mock

import server
from snaclex import uniprot


def _fake_record(acc, gene, taxon, organism, length):
    return {
        "canonical_id": f"SNX:PRT:{taxon}:{acc}:1",
        "accessions": {"uniprot_primary": acc},
        "preferred_name": f"{gene} protein", "gene": {"symbol": gene},
        "organism": {"scientific_name": organism, "taxon_id": taxon},
        "review_status": "reviewed", "sequence": {"length": length},
        "sequence_analysis": {"length": length, "molecular_weight_Da": length * 110.0,
                              "isoelectric_point": 6.5, "gravy": -0.2},
        "domains_motifs": [{"name": "Protein kinase"}], "variants": [],
        "structures": {"tier": "experimental"},
    }

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

    def _post(self, path, payload):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("POST", path, body=json.dumps(payload),
                     headers={"Content-Type": "application/json"})
        resp = conn.getresponse()
        body = resp.read()
        conn.close()
        return resp, json.loads(body)

    def _poll(self, job_id, tries=100):
        for _ in range(tries):
            resp, out = self._get(f"/api/jobs/{job_id}")
            if out.get("status") in ("done", "error"):
                return out
            time.sleep(0.02)
        return out

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

    def test_run_batch_job_direct(self):
        recs = {
            "P00519": _fake_record("P00519", "ABL1", 9606, "Homo sapiens", 1130),
            "P00520": _fake_record("P00520", "ABL1", 10090, "Mus musculus", 1123),
        }

        def fake_build(q, **kw):
            return (recs[q], None) if q in recs else (None, "not found")

        with mock.patch.object(server, "build_protein_record", side_effect=fake_build):
            out = server.run_batch_job(
                {"queries": ["P00519", "P00520", "NOPE"], "compare": True})
        self.assertEqual(out["n_requested"], 3)
        self.assertEqual(out["n_resolved"], 2)
        self.assertEqual(len(out["errors"]), 1)
        self.assertIn("ABL1", out["comparison"]["orthology"]["ortholog_candidates"])

    def test_batch_endpoint_submits_and_completes(self):
        rec = _fake_record("P00519", "ABL1", 9606, "Homo sapiens", 1130)
        with mock.patch.object(server, "build_protein_record",
                               return_value=(rec, None)):
            resp, out = self._post("/api/protein/batch",
                                   {"queries": ["P00519"], "compare": True})
            self.assertEqual(resp.status, 202)
            job_id = out["job_id"]
            result = self._poll(job_id)
        self.assertEqual(result["status"], "done")
        self.assertEqual(result["result"]["n_resolved"], 1)
        self.assertIn("comparison", result["result"])

    def test_batch_endpoint_rejects_empty(self):
        resp, out = self._post("/api/protein/batch", {"queries": []})
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

    def test_homologs_endpoint(self):
        res = _resolution_with_uniprot()
        res["_uniprot_fragment"]["sequence"]["value"] = "A" * 40
        homo = {"available": True, "method": {"identity_cutoff": 0.3},
                "structures": [{"pdb_id": "2GQG", "sequence_identity": 0.62,
                                "docking_suitable": True}]}
        with mock.patch.object(server.idresolve, "resolve", return_value=res), \
             mock.patch.object(server.homology, "find_homologous_structures",
                               return_value=homo):
            resp, out = self._get("/api/homologs?acc=P00519")
        self.assertEqual(resp.status, 200)
        self.assertTrue(out["available"])
        self.assertEqual(out["structures"][0]["pdb_id"], "2GQG")

    def test_variant_maps_residue(self):
        res = _resolution_with_uniprot()
        res["_uniprot_fragment"]["sequence"]["value"] = "MDLSAKLICE"
        res["_uniprot_fragment"]["domains_motifs"] = [
            {"name": "Kinase", "type": "domain", "start": 1, "end": 10}]
        avail = {"experimental": [{"pdb_id": "1ABC"}], "homologous": [],
                 "predicted": [], "tier": "experimental",
                 "best_for_docking": "1ABC", "warnings": []}
        with mock.patch.object(server.idresolve, "resolve", return_value=res), \
             mock.patch.object(server.proteinrecord, "structure_availability",
                               return_value=avail):
            resp, out = self._get("/api/variant?q=GENE+K6R")
        self.assertEqual(resp.status, 200)
        self.assertEqual(out["variant"]["notation"], "K6R")
        self.assertTrue(out["analysis"]["sequence_mapping"]["wt_matches"])
        self.assertTrue(out["analysis"]["domain_disruption"])
        self.assertEqual(out["analysis"]["consequence"], "missense")
        self.assertIn("Not a clinical interpretation",
                      out["analysis"]["clinical"]["disclaimer"])
        self.assertEqual(out["analysis"]["structure_mapping"]["experimental_structures"],
                         ["1ABC"])

    def test_variant_bad_input_400(self):
        resp, out = self._get("/api/variant?q=not+a+variant")
        self.assertEqual(resp.status, 400)

    def test_workspace_oncology_lens(self):
        res = _resolution_with_uniprot()
        res["_uniprot_fragment"]["keywords"] = [
            {"name": "Proto-oncogene", "category": "x"}]
        with mock.patch.object(server.idresolve, "resolve", return_value=res):
            resp, out = self._get("/api/workspace?acc=P00519&view=oncology")
        self.assertEqual(resp.status, 200)
        self.assertIn("oncology", out)
        self.assertNotIn("immunology", out)
        roles = {r["role"] for r in out["oncology"]["roles"]}
        self.assertIn("oncogene", roles)
        self.assertIn("NO claim", out["oncology"]["separation_note"])

    def test_workspace_bad_view_400(self):
        resp, out = self._get("/api/workspace?acc=P00519&view=nonsense")
        # resolve isn't reached; view validation fails first only if resolve ok —
        # but bad view is rejected before building, so expect 400.
        self.assertIn(resp.status, (400,))

    def test_evidence_level_d_candidates_opt_in(self):
        res = _resolution_with_uniprot()
        res["_uniprot_fragment"]["sequence"]["value"] = "A" * 40
        homo = {"available": True, "method": {},
                "structures": [{"pdb_id": "2GQG", "sequence_identity": 0.62}]}
        with mock.patch.object(server.idresolve, "resolve", return_value=res), \
             mock.patch.object(server.pubchem, "lookup_compound", return_value={"cid": 5291}), \
             mock.patch.object(server.pubchem, "fetch_identity",
                               return_value={"cid": 5291, "inchikey": "X"}), \
             mock.patch.object(server.pubchem, "fetch_parent_cid", return_value=None), \
             mock.patch.object(server.pubchem, "bioassay_summary", return_value=[]), \
             mock.patch.object(server.proteinrecord, "structure_availability",
                               return_value={"tier": "homologous", "best_for_docking": "2GQG"}), \
             mock.patch.object(server.homology, "find_homologous_structures",
                               return_value=homo):
            resp, out = self._get("/api/evidence?protein=P00519&chemical=imatinib&homologs=1")
        self.assertEqual(resp.status, 200)
        self.assertIsNotNone(out["level_d_candidates"])
        self.assertEqual(out["level_d_candidates"]["structural_homologs"][0]["pdb_id"], "2GQG")
        # A Level-D "binds" claim is NOT asserted — only candidates are surfaced.
        self.assertNotIn("D", out["levels_not_gathered_in_slice"])


if __name__ == "__main__":
    unittest.main()
