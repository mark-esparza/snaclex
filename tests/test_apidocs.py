"""Tests for the API contract served at /api/docs."""

import json
import unittest

from snaclex import __version__, apidocs


class TestApiDocs(unittest.TestCase):
    def setUp(self):
        self.c = apidocs.contract()

    def test_top_level_shape(self):
        for key in ("tool", "version", "limits", "errors", "endpoints"):
            self.assertIn(key, self.c)
        self.assertEqual(self.c["version"], __version__)

    def test_every_endpoint_well_formed(self):
        for ep in self.c["endpoints"]:
            self.assertIn(ep["method"], ("GET", "POST"))
            self.assertTrue(ep["path"].startswith("/api/"))
            self.assertIn("returns", ep)

    def test_key_endpoints_present(self):
        paths = {e["path"] for e in self.c["endpoints"]}
        for p in ("/api/analyze", "/api/jobs", "/api/jobs/{id}",
                  "/api/upload", "/api/docs", "/api/version"):
            self.assertIn(p, paths)

    def test_serializable(self):
        json.dumps(self.c)


if __name__ == "__main__":
    unittest.main()
