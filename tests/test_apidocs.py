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

    def test_every_routed_path_is_documented(self):
        """Guard against the contract drifting behind the router.

        Scrapes the literal paths the GET router compares against and asserts
        each has a contract entry, so a new endpoint cannot ship undocumented.
        """
        import os
        import re

        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        with open(os.path.join(root, "server.py"), encoding="utf-8") as fh:
            source = fh.read()

        routed = set(re.findall(r'path == "(/api/[^"]*)"', source))
        routed |= {
            prefix + "{id}"
            for prefix in re.findall(r'path\.startswith\("(/api/[^"]*/)"\)', source)
        }
        documented = {e["path"] for e in self.c["endpoints"]}

        missing = {
            path for path in routed
            if path not in documented
            and path.replace("{id}", "") not in
            {d.split("{")[0] for d in documented}
        }
        self.assertEqual(missing, set(), f"undocumented endpoints: {missing}")


if __name__ == "__main__":
    unittest.main()
