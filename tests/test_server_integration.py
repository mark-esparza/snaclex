"""End-to-end tests against a live in-process server.

Boots the real ThreadingHTTPServer on an ephemeral port and exercises endpoints
that need no upstream network (version, static, error paths) plus the security
header pipeline.
"""

import http.client
import json
import threading
import time
import unittest
from http.server import ThreadingHTTPServer
from unittest import mock

import server


class TestServerIntegration(unittest.TestCase):
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
        # Give each test a fresh rate-limit budget (shared module singletons).
        server._IP_LIMITER._buckets.clear()
        server._HEAVY_LIMITER._buckets.clear()

    def _get(self, path, headers=None):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("GET", path, headers=headers or {})
        resp = conn.getresponse()
        body = resp.read()
        conn.close()
        return resp, body

    def _post(self, path, payload):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        body = json.dumps(payload)
        conn.request("POST", path, body=body,
                     headers={"Content-Type": "application/json"})
        resp = conn.getresponse()
        data = resp.read()
        conn.close()
        return resp, data

    def test_version_endpoint(self):
        resp, body = self._get("/api/version")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertEqual(data["version"], server.SNACLEX_VERSION)
        self.assertTrue(data["research_only"])

    def test_security_headers_on_static(self):
        resp, _ = self._get("/")
        self.assertEqual(resp.status, 200)
        self.assertEqual(resp.getheader("X-Content-Type-Options"), "nosniff")
        self.assertEqual(resp.getheader("X-Frame-Options"), "DENY")
        self.assertIn("3Dmol.org", resp.getheader("Content-Security-Policy"))

    def test_hsts_only_behind_https_proxy(self):
        resp, _ = self._get("/api/version")
        self.assertIsNone(resp.getheader("Strict-Transport-Security"))
        resp2, _ = self._get("/api/version", {"X-Forwarded-Proto": "https"})
        self.assertIsNotNone(resp2.getheader("Strict-Transport-Security"))

    def test_static_scope_blocks_non_web_files(self):
        # server.py lives at the repo root, not under web/ — must not be served.
        resp, _ = self._get("/server.py")
        self.assertEqual(resp.status, 404)

    def test_missing_param_returns_400_json(self):
        resp, body = self._get("/api/analyze")
        self.assertEqual(resp.status, 400)
        self.assertIn("error", json.loads(body))

    # ---- async job queue (stubbed runner, no upstream network) -------
    def test_job_lifecycle_done(self):
        with mock.patch.dict(server._JOB_RUNNERS,
                             {"echo": lambda p: {"got": p}}, clear=False):
            resp, body = self._post("/api/jobs",
                                    {"kind": "echo", "params": {"hi": 1}})
            self.assertEqual(resp.status, 202)
            job_id = json.loads(body)["job_id"]

            for _ in range(50):
                r, b = self._get(f"/api/jobs/{job_id}")
                st = json.loads(b)
                if st["status"] == "done":
                    self.assertEqual(st["result"], {"got": {"hi": 1}})
                    break
                time.sleep(0.02)
            else:
                self.fail("job never completed")

    def test_job_error_is_reported(self):
        def boom(_p):
            raise ValueError("nope")

        with mock.patch.dict(server._JOB_RUNNERS, {"boom": boom}, clear=False):
            resp, body = self._post("/api/jobs", {"kind": "boom", "params": {}})
            job_id = json.loads(body)["job_id"]
            for _ in range(50):
                _r, b = self._get(f"/api/jobs/{job_id}")
                st = json.loads(b)
                if st["status"] == "error":
                    self.assertIn("nope", st["error"])
                    break
                time.sleep(0.02)
            else:
                self.fail("job error never surfaced")

    def test_unknown_job_kind_400(self):
        resp, body = self._post("/api/jobs", {"kind": "nonsense", "params": {}})
        self.assertEqual(resp.status, 400)
        self.assertIn("error", json.loads(body))

    def test_unknown_job_id_404(self):
        resp, _ = self._get("/api/jobs/deadbeef")
        self.assertEqual(resp.status, 404)

    # ---- API docs ----------------------------------------------------
    def test_api_docs(self):
        resp, body = self._get("/api/docs")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertIn("endpoints", data)
        self.assertTrue(any(e["path"] == "/api/upload" for e in data["endpoints"]))

    # ---- structure upload --------------------------------------------
    def _post_raw(self, path, text, ctype="text/plain"):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("POST", path, body=text, headers={"Content-Type": ctype})
        resp = conn.getresponse()
        data = resp.read()
        conn.close()
        return resp, data

    def test_upload_and_analyze(self):
        resp, body = self._post_raw("/api/upload", _SAMPLE_PDB)
        self.assertEqual(resp.status, 200)
        up = json.loads(body)
        self.assertTrue(up["upload_id"].startswith("UL"))
        self.assertGreaterEqual(up["protein_atom_count"], 10)
        # The uploaded structure is now loadable by its id like any PDB.
        resp2, body2 = self._get(f"/api/analyze?pdb={up['upload_id']}")
        self.assertEqual(resp2.status, 200)
        self.assertEqual(json.loads(body2)["protein_atom_count"],
                         up["protein_atom_count"])

    def test_upload_rejects_junk(self):
        resp, body = self._post_raw("/api/upload", "not a structure at all\n")
        self.assertEqual(resp.status, 400)
        self.assertIn("error", json.loads(body))

    def test_upload_accepts_mmcif(self):
        from tests.fixtures import mmcif_text
        rows = [{"element": "C", "name": f"C{i}", "seq": i, "x": float(i)}
                for i in range(12)]
        body_text = mmcif_text(rows)
        resp, body = self._post_raw("/api/upload", body_text)
        self.assertEqual(resp.status, 200)
        up = json.loads(body)
        self.assertTrue(up["upload_id"].startswith("UL"))
        # The viewer payload is re-serialized to PDB format, not raw mmCIF.
        self.assertNotIn("_atom_site.", up["pdb_data"])
        self.assertIn("ATOM", up["pdb_data"])

    def test_upload_rejects_oversized(self):
        from tests.fixtures import mmcif_text
        rows = [{"element": "C", "name": f"C{i}", "seq": i} for i in range(12)]
        with mock.patch.object(server, "MAX_STRUCTURE_ATOMS", 5):
            resp, body = self._post_raw("/api/upload", mmcif_text(rows))
        self.assertEqual(resp.status, 413)
        self.assertIn("too large", json.loads(body)["error"])

    # ---- large assembly -> per-chain loading -------------------------
    def _two_chain_structure(self):
        from tests.fixtures import atom, structure
        prot = (
            [atom("C", float(i), 0, 0, name=f"A{i}", res_name="LEU", chain="A", res_seq=i)
             for i in range(6)]
            + [atom("C", float(i), 5, 0, name=f"B{i}", res_name="LEU", chain="B", res_seq=i)
               for i in range(6)]
        )
        return structure(prot)

    def test_analyze_too_large_offers_chains(self):
        s = self._two_chain_structure()
        meta = {"pdb_id": "TST1", "title": "Test assembly"}
        with mock.patch.object(server, "_load_full", return_value=("ATOMS", s, meta)), \
             mock.patch.object(server, "MAX_STRUCTURE_ATOMS", 8):
            resp, body = self._get("/api/analyze?pdb=TST1")
        data = json.loads(body)
        self.assertEqual(resp.status, 200)
        self.assertTrue(data["too_large"])
        self.assertEqual(data["n_atoms"], 12)
        self.assertEqual({c["chain"] for c in data["chains"]}, {"A", "B"})

    def test_analyze_single_chain_loads(self):
        s = self._two_chain_structure()
        meta = {"pdb_id": "TST2", "title": "Test assembly"}
        with mock.patch.object(server, "_load_full", return_value=("ATOMS", s, meta)), \
             mock.patch.object(server, "MAX_STRUCTURE_ATOMS", 8):
            resp, body = self._get("/api/analyze?pdb=TST2&chain=A")
        data = json.loads(body)
        self.assertEqual(resp.status, 200)
        self.assertEqual(data["id"], "TST2-A")          # synthetic subset id
        self.assertEqual(data["protein_atom_count"], 6)  # only chain A
        self.assertNotIn("_atom_site.", data["pdb_data"])
        # The subset is now loadable by its id like any structure.
        resp2, _ = self._get("/api/analyze?pdb=TST2-A")
        self.assertEqual(resp2.status, 200)

    def test_analyze_unknown_chain_404(self):
        s = self._two_chain_structure()
        meta = {"pdb_id": "TST3", "title": "x"}
        with mock.patch.object(server, "_load_full", return_value=("ATOMS", s, meta)):
            resp, _ = self._get("/api/analyze?pdb=TST3&chain=Z")
        self.assertEqual(resp.status, 404)

    # ---- Benchmark Mode ----------------------------------------------
    def test_benchmark_cases_listed(self):
        resp, body = self._get("/api/benchmark/cases")
        self.assertEqual(resp.status, 200)
        cases = json.loads(body)["cases"]
        self.assertTrue(cases)
        self.assertTrue(all("pdb" in c for c in cases))

    def _structure_with_ligand(self):
        from tests.fixtures import atom, component, structure
        prot = [atom("N", 3, 0, 0, name="N", res_name="ALA", chain="A", res_seq=i)
                for i in range(12)]
        lig = component("LIG", [
            atom("O", 0, 0, 0, hetero=True, res_name="LIG", chain="B", res_seq=900),
            atom("C", 1.4, 0, 0, hetero=True, res_name="LIG", chain="B", res_seq=900),
            atom("C", 1.4, 1.4, 0, hetero=True, res_name="LIG", chain="B", res_seq=900),
        ], chain="B", res_seq=900)
        return structure(prot, [lig])

    def test_benchmark_job(self):
        s = self._structure_with_ligand()
        meta = {"pdb_id": "BENCH", "title": "Bench case"}
        with mock.patch.object(server, "_load_structure",
                               return_value=("ATOMS", s, meta)):
            resp, body = self._post("/api/jobs", {"kind": "benchmark",
                                                  "params": {"pdb": "BENCH"}})
            self.assertEqual(resp.status, 202)
            job_id = json.loads(body)["job_id"]
            for _ in range(200):
                _r, b = self._get(f"/api/jobs/{job_id}")
                st = json.loads(b)
                if st["status"] == "done":
                    res = st["result"]
                    self.assertIn("pocket", res)
                    self.assertIn("physical_plausibility", res)
                    self.assertIn("interactions_total", res)
                    break
                if st["status"] == "error":
                    self.fail(f"benchmark job errored: {st['error']}")
                time.sleep(0.05)
            else:
                self.fail("benchmark job did not finish")


    # ---- interface / nucleic / HLA / esm (mocked structures) ---------
    def _interface_structure(self):
        from tests.fixtures import atom, structure
        prot = [
            atom("N", 0, 0, 0, name="NH1", res_name="ARG", chain="A", res_seq=10),
            atom("O", 3.0, 0, 0, name="OD1", res_name="ASP", chain="B", res_seq=20),
        ]
        return structure(prot)

    def test_interface_endpoint(self):
        s = self._interface_structure()
        meta = {"pdb_id": "IFC1", "title": "complex"}
        with mock.patch.object(server, "_load_structure",
                               return_value=("ATOMS", s, meta)):
            resp, body = self._get("/api/interface?pdb=IFC1&a=A&b=B")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertEqual(data["profile"]["mode"], "protein-protein")
        self.assertGreaterEqual(data["profile"]["counts"]["salt_bridge"], 1)
        self.assertIn("methods", data)

    def test_nucleic_endpoint(self):
        from tests.fixtures import atom, structure
        prot = [atom("N", 0, 0, 0, name="NZ", res_name="LYS", chain="A", res_seq=10)]
        nuc = [atom("O", 3.0, 0, 0, name="OP1", res_name="DA", chain="B", res_seq=1)]
        s = structure(prot, nucleic=nuc)
        meta = {"pdb_id": "NUC1", "title": "protein-DNA"}
        with mock.patch.object(server, "_load_structure",
                               return_value=("ATOMS", s, meta)):
            resp, body = self._get("/api/nucleic?pdb=NUC1")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertTrue(data["available"])
        self.assertEqual(data["profile"]["mode"], "protein-nucleic")

    def test_nucleic_endpoint_no_nucleic(self):
        s = self._interface_structure()
        meta = {"pdb_id": "NUC2", "title": "no DNA"}
        with mock.patch.object(server, "_load_structure",
                               return_value=("ATOMS", s, meta)):
            resp, body = self._get("/api/nucleic?pdb=NUC2")
        self.assertEqual(resp.status, 200)
        self.assertFalse(json.loads(body)["available"])

    def _hla_structure(self):
        from tests.fixtures import atom, structure
        groove = {100, 150}
        prot = []
        for i in range(1, 271):
            if i in groove:
                continue
            prot.append(atom("C", 500.0 + i, 0, 0, name="CA", res_name="ALA",
                             chain="A", res_seq=i))
        prot.append(atom("N", 3.0, 0, 0, name="N", res_name="GLN", chain="A", res_seq=100))
        prot.append(atom("O", 3.0, 20, 0, name="OD1", res_name="ASP", chain="A", res_seq=150))
        for i in range(1, 100):
            prot.append(atom("C", -500.0 - i, 0, 0, name="CA", res_name="ALA",
                             chain="B", res_seq=i))
        pep = {2: (0.0, 0.0, 0.0), 9: (0.0, 20.0, 0.0)}
        for i in range(1, 10):
            x, y, z = pep.get(i, (100.0, 100.0 + i, 0.0))
            prot.append(atom("O", x, y, z, name="O", res_name="GLY", chain="C", res_seq=i))
        return structure(prot)

    def test_hla_endpoint(self):
        s = self._hla_structure()
        meta = {"pdb_id": "HLA1", "title": "HLA-A*02:01 complex"}
        with mock.patch.object(server, "_load_structure",
                               return_value=("ATOMS", s, meta)), \
             mock.patch.object(server, "_get_uniprots", return_value=[]):
            resp, body = self._get("/api/hla?pdb=HLA1")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertTrue(data["detection"]["is_hla"])
        self.assertEqual(data["detection"]["mhc_class"], "I")
        self.assertTrue(data["groove"]["available"])
        self.assertEqual(data["groove"]["anchor_peptide_positions"], [2, 9])

    def test_hla_cases_endpoint(self):
        resp, body = self._get("/api/hla/cases")
        self.assertEqual(resp.status, 200)
        cases = json.loads(body)["cases"]
        self.assertTrue(any(c["allele"].startswith("HLA-B*57") for c in cases))

    def test_esm_status_endpoint(self):
        resp, body = self._get("/api/esm")
        self.assertEqual(resp.status, 200)
        data = json.loads(body)
        self.assertIn("available", data)
        self.assertIn("fold_model", data)

    def _seq_structure(self):
        from tests.fixtures import atom, structure
        three = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
                 "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU"}
        seq = "ACDEFGHIKL"
        prot = [atom("C", float(i), 0, 0, name="CA", res_name=three[a],
                     chain="A", res_seq=i + 1) for i, a in enumerate(seq)]
        return structure(prot)

    def test_variants_job(self):
        s = self._seq_structure()
        meta = {"pdb_id": "VAR1", "title": "p53-like"}
        with mock.patch.object(server, "_load_structure",
                               return_value=("ATOMS", s, meta)), \
             mock.patch.object(server, "_get_uniprots", return_value=["P99999"]), \
             mock.patch.object(server.variants, "fetch_uniprot_sequence",
                               return_value="ACDEFGHIKL"):
            resp, body = self._post("/api/jobs", {"kind": "variants",
                                                  "params": {"pdb": "VAR1",
                                                             "variants": ["C2D", "K9A"]}})
            self.assertEqual(resp.status, 202)
            job_id = json.loads(body)["job_id"]
            for _ in range(200):
                _r, b = self._get(f"/api/jobs/{job_id}")
                st = json.loads(b)
                if st["status"] == "done":
                    res = st["result"]
                    self.assertEqual(res["mapping"]["mapped_count"], 2)
                    self.assertFalse(res["esm"]["available"])
                    break
                if st["status"] == "error":
                    self.fail(f"variants job errored: {st['error']}")
                time.sleep(0.05)
            else:
                self.fail("variants job did not finish")


    # ---- component chemical-name enrichment --------------------------
    def _ligand_structure(self, code):
        from tests.fixtures import atom, component, structure
        lig = component(code, [
            atom("C", 0, 0, 0, hetero=True, res_name=code, chain="B", res_seq=900),
            atom("N", 1.4, 0, 0, hetero=True, res_name=code, chain="B", res_seq=900),
        ], chain="B", res_seq=900)
        prot = [atom("C", 9, 0, 0, name="CA", res_name="ALA", chain="A", res_seq=1)]
        return structure(prot, [lig])

    def test_components_enriched_with_chem_name(self):
        s = self._ligand_structure("STI")
        server._CHEMCOMP_CACHE.clear()
        ccd = {"STI": {"name": "Imatinib", "formula": "C29 H31 N7 O",
                       "formula_weight": 493.6, "smiles": "Cc1...", "synonyms": ["Gleevec"]}}
        with mock.patch.object(server.rcsb, "fetch_chem_components", return_value=ccd):
            comps = server._components_json(s)
        sti = next(c for c in comps if c["res_name"] == "STI")
        self.assertEqual(sti["chem_name"], "Imatinib")
        self.assertEqual(sti["formula"], "C29 H31 N7 O")
        self.assertEqual(sti["smiles"], "Cc1...")

    def test_components_graceful_without_ccd(self):
        s = self._ligand_structure("XYZ")
        server._CHEMCOMP_CACHE.clear()
        with mock.patch.object(server.rcsb, "fetch_chem_components", return_value={}):
            comps = server._components_json(s)
        xyz = next(c for c in comps if c["res_name"] == "XYZ")
        self.assertIsNone(xyz["chem_name"])
        self.assertEqual(xyz["res_name"], "XYZ")  # still present, just unnamed


def _pdb_line(rec, serial, name, res, chain, seq, x, y, z, el):
    return (f"{rec:<6}{serial:>5} {name:<4} {res:>3} {chain}{seq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}")


# A minimal but column-correct PDB with 12 protein atoms.
_SAMPLE_PDB = "\n".join(
    _pdb_line("ATOM", i + 1, el, "LEU", "A", i + 1, float(i), 0.0, 0.0, el)
    for i, el in enumerate(["N", "C", "C", "O", "C", "C",
                            "N", "C", "C", "O", "C", "C"])
) + "\nEND\n"


if __name__ == "__main__":
    unittest.main()
