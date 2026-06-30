"""Headless browser smoke test of the golden analysis path.

This is the verification the follow-up audit flagged as P0: ~2,400 lines of
front-end JavaScript that, until now, had only ever been syntax-checked. It
boots the real server and drives a real Chromium through:

    load a structure -> pick a bound molecule -> see interactions

while asserting **no uncaught JS error** fires — the single highest-value check
for catching render-path regressions.

It runs fully offline: the upstream RCSB calls are served from a pre-seeded disk
cache (no network), so CI is deterministic. Requires Playwright + Chromium
(see requirements-dev.txt); skipped automatically if Playwright isn't installed.
"""

from __future__ import annotations

import os
import socket
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.request

try:
    from playwright.sync_api import expect, sync_playwright
    _HAVE_PLAYWRIGHT = True
except ImportError:  # keep importable (and skippable) without the dev dep
    _HAVE_PLAYWRIGHT = False

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_ID = "1ABC"  # served from the seeded cache, not the real PDB


def _pdb_line(rec, serial, name, res, chain, seq, x, y, z, el):
    return (f"{rec:<6}{serial:>5} {name:<4} {res:>3} {chain}{seq:>4}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}")


def _fixture_pdb():
    rows = []
    serial = 1
    # 12 protein atoms across three residues (chain A).
    for r, base_y in ((10, 0.0), (11, 3.0), (12, -3.0)):
        for nm, el, dx in (("N", "N", 3.0), ("CA", "C", 4.0),
                           ("C", "C", 5.0), ("O", "O", 6.0)):
            rows.append(_pdb_line("ATOM", serial, nm, "ALA", "A", r, dx, base_y, 0.0, el))
            serial += 1
    # A 3-atom ligand whose O sits 3.0 A from the res-10 backbone N -> H-bond.
    for nm, el, (x, y, z) in (("O1", "O", (0.0, 0.0, 0.0)),
                              ("C1", "C", (1.4, 0.0, 0.0)),
                              ("C2", "C", (1.4, 1.4, 0.0))):
        rows.append(_pdb_line("HETATM", serial, nm, "LIG", "B", 900, x, y, z, el))
        serial += 1
    rows.append("END")
    return "\n".join(rows) + "\n"


def _seed_cache(directory):
    """Pre-seed the disk HTTP cache so the server needs no network."""
    sys.path.insert(0, REPO)
    from snaclex.cache import DiskCache
    cache = DiskCache(directory, ttl_seconds=10 * 365 * 86400)
    cache.set(f"https://files.rcsb.org/download/{TEST_ID}.pdb",
              _fixture_pdb().encode())
    cache.set(f"https://data.rcsb.org/rest/v1/core/entry/{TEST_ID}",
              b'{"struct": {"title": "E2E test structure"}}')


def _free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


@unittest.skipUnless(_HAVE_PLAYWRIGHT, "playwright not installed")
class TestGoldenPath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache_dir = tempfile.mkdtemp(prefix="snaclex-e2e-cache-")
        _seed_cache(cls.cache_dir)
        cls.port = _free_port()
        env = dict(os.environ, SNACLEX_HTTP_CACHE=cls.cache_dir)
        cls.proc = subprocess.Popen(
            [sys.executable, "server.py", "--port", str(cls.port), "--host", "127.0.0.1"],
            cwd=REPO, env=env,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        cls._wait_until_up()
        cls.pw = sync_playwright().start()
        cls.browser = cls.pw.chromium.launch()

    @classmethod
    def tearDownClass(cls):
        try:
            cls.browser.close()
            cls.pw.stop()
        finally:
            cls.proc.terminate()
            try:
                cls.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                cls.proc.kill()

    @classmethod
    def _wait_until_up(cls, timeout=20):
        url = f"http://127.0.0.1:{cls.port}/api/version"
        end = time.time() + timeout
        while time.time() < end:
            try:
                with urllib.request.urlopen(url, timeout=2):
                    return
            except OSError:
                time.sleep(0.2)
        raise RuntimeError("server did not start")

    def test_load_pick_and_profile_without_js_errors(self):
        page = self.browser.new_page()
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))

        page.goto(f"http://127.0.0.1:{self.port}/")
        self.assertIn("SnaCleX", page.title())

        page.fill("#pdbInput", TEST_ID)
        page.click("#loadBtn")
        # Use expect() assertions (which poll over the CDP protocol) rather than
        # wait_for_function — the latter evals a JS string in the page, which the
        # production Content-Security-Policy correctly blocks (no 'unsafe-eval').
        expect(page.locator("#statusBar")).to_contain_text("Loaded", timeout=15000)

        # A bound molecule (the ligand) should be listed; pick it.
        page.wait_for_selector(".comp", timeout=5000)
        page.locator(".comp").first.click()

        # The interactions panel should leave its empty placeholder.
        expect(page.locator("#interactionContent")).not_to_contain_text(
            "Select a bound molecule", timeout=15000
        )

        # The 3D viewer's self-describing legend renders on load and explains the
        # visual encoding (verifies renderLegend runs without a JS error).
        page.click('.tab[data-tab="viewer"]')
        expect(page.locator("#viewerLegend")).to_contain_text("Protein", timeout=10000)

        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()

    def test_new_tabs_render_without_js_errors(self):
        # Drives the Interface / Variants / HLA render paths added for the
        # nucleic-acid, antibody/antigen and HLA features. The seeded fixture is
        # a single-chain protein, so these exercise the no-result / not-detected
        # branches — the point is to catch JS render-path errors in the new code.
        page = self.browser.new_page()
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))

        page.goto(f"http://127.0.0.1:{self.port}/")
        page.fill("#pdbInput", TEST_ID)
        page.click("#loadBtn")
        expect(page.locator("#statusBar")).to_contain_text("Loaded", timeout=15000)

        page.click('.tab[data-tab="interface"]')
        page.click("#interfaceBtn")

        page.click('.tab[data-tab="hla"]')
        page.click("#hlaBtn")
        expect(page.locator("#hlaContent")).to_contain_text("detected", timeout=15000)

        page.click('.tab[data-tab="variants"]')  # just switch + render placeholder

        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()


if __name__ == "__main__":
    unittest.main()
