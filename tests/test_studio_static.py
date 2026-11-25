"""The Animation Studio page: served, linked, bare, and dependency-safe.

The studio is static files under ``web/``. These tests pin the parts the
server and the rest of the site are responsible for: the page and its example
scripts are served with the right types, the vendored libraries never import a
bare module name (the CSP forbids the import map that would resolve one), every
page links to the studio, and the editor opens empty (the same no-suggestions
rule as the workbench).
"""

import http.client
import os
import re
import threading
import unittest
from http.server import ThreadingHTTPServer

import server

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WEB = os.path.join(ROOT, "web")
PAGES = ("index.html", "help.html", "api.html", "privacy.html", "terms.html")
DEMO_PDB_IDS = ("1HSG", "1CA2", "4HHB", "3PTB", "4DFR", "6BCX")


def _read(*parts):
    with open(os.path.join(WEB, *parts), encoding="utf-8") as fh:
        return fh.read()


class TestStudioServed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
        cls.port = cls.httpd.server_address[1]
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.httpd.shutdown()

    def setUp(self):
        server._IP_LIMITER._buckets.clear()

    def _get(self, path):
        conn = http.client.HTTPConnection("127.0.0.1", self.port, timeout=5)
        conn.request("GET", path)
        resp = conn.getresponse()
        body = resp.read()
        conn.close()
        return resp, body

    def test_studio_page_served(self):
        resp, body = self._get("/studio.html")
        self.assertEqual(resp.status, 200)
        self.assertTrue(resp.getheader("Content-Type").startswith("text/html"))
        text = body.decode("utf-8")
        self.assertIn('type="module"', text)
        self.assertIn("studio/app.js", text)

    def test_example_txt_served_as_text(self):
        resp, _ = self._get("/studio/examples/transcription-basic.txt")
        self.assertEqual(resp.status, 200)
        self.assertTrue(resp.getheader("Content-Type").startswith("text/plain"))


class TestStudioFiles(unittest.TestCase):
    def test_vendor_three_has_no_bare_imports(self):
        vendor = os.path.join(WEB, "vendor")
        files = [f for f in os.listdir(vendor) if f.endswith(".js")]
        self.assertTrue(files, "web/vendor has no .js files")
        for name in files:
            text = _read("vendor", name)
            self.assertNotRegex(text, r"""from\s*['"]three['"]""", name)

    def test_studio_linked_from_every_page(self):
        for name in PAGES:
            self.assertIn('href="/studio.html"', _read(name), name)

    def test_editor_starts_empty(self):
        html = _read("studio.html")
        m = re.search(r"<textarea[^>]*id=\"editor\"[^>]*>(.*?)</textarea>", html, re.S)
        self.assertIsNotNone(m, "no <textarea id=\"editor\">")
        self.assertEqual(m.group(1).strip(), "")
        for pid in DEMO_PDB_IDS:
            self.assertNotIn(pid, html)


if __name__ == "__main__":
    unittest.main()
