"""Headless browser smoke tests of the Animation Studio page.

Boots the real server (offline: RCSB responses come from a pre-seeded disk
cache, shared with test_smoke) and drives Chromium through the studio:
loading an example, playing it, recovering from a half-typed script, seeking,
and building a protein from a real (fixture) structure — always asserting that
no uncaught JS error fires.

The page exposes ``window.__studio`` only when loaded with ``?test=1``.
Polling is done from Python with page.evaluate(): the production CSP blocks
the string-eval that page.wait_for_function() relies on.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import time
import unittest

from e2e.test_smoke import _HAVE_PLAYWRIGHT, REPO, TEST_ID, _free_port, _seed_cache

if _HAVE_PLAYWRIGHT:
    from playwright.sync_api import expect, sync_playwright


def _poll(page, js, timeout=15.0, interval=0.1):
    end = time.time() + timeout
    value = None
    while time.time() < end:
        value = page.evaluate(js)
        if value:
            return value
        time.sleep(interval)
    return value


@unittest.skipUnless(_HAVE_PLAYWRIGHT, "playwright not installed")
class TestStudio(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache_dir = tempfile.mkdtemp(prefix="snaclex-e2e-studio-")
        _seed_cache(cls.cache_dir)
        cls.port = _free_port()
        env = dict(os.environ, SNACLEX_HTTP_CACHE=cls.cache_dir)
        cls.proc = subprocess.Popen(
            [sys.executable, "server.py", "--port", str(cls.port), "--host", "127.0.0.1"],
            cwd=REPO, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        import urllib.request
        end = time.time() + 20
        while True:
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{cls.port}/api/version", timeout=2):
                    break
            except OSError:
                if time.time() > end:
                    raise RuntimeError("server did not start")
                time.sleep(0.2)
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

    def _open(self):
        page = self.browser.new_page(viewport={"width": 1280, "height": 900})
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(f"http://127.0.0.1:{self.port}/studio.html?test=1")
        self.assertTrue(_poll(page, "() => !!window.__studio"), "studio did not start")
        return page, errors

    def _load_example(self, page, value="transcription-basic"):
        page.select_option("#examples", value)
        self.assertTrue(
            _poll(page, "() => { const t = window.__studio.getTimeline(); return !!t && t.order.length > 0; }"),
            "example did not compile",
        )

    def test_example_loads_and_plays(self):
        page, errors = self._open()
        self._load_example(page)
        self.assertNotIn("Line", page.locator("#diagnostics .err").all_inner_texts().__str__())
        page.click("#play")
        page.wait_for_timeout(1000)
        self.assertTrue(page.evaluate("() => window.__studio.isPlaying()"))
        self.assertGreater(page.evaluate("() => window.__studio.time()"), 0.3)
        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()

    def test_broken_script_keeps_last_good_preview(self):
        page, errors = self._open()
        self._load_example(page)
        names = page.evaluate("() => window.__studio.getTimeline().order")
        page.locator("#editor").press("Control+End")
        page.locator("#editor").type("\n3s dna unwnd")
        expect(page.locator("#diagnostics")).to_contain_text("did you mean", timeout=5000)
        self.assertEqual(page.evaluate("() => window.__studio.getTimeline().order"), names)
        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()

    def test_seek_updates_clock(self):
        page, errors = self._open()
        self._load_example(page)
        page.evaluate("() => window.__studio.seek(2)")
        expect(page.locator("#clock")).to_contain_text("0:02.0", timeout=3000)
        self.assertTrue(page.locator("#clock").inner_text().startswith("0:02.0"))
        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()

    def test_protein_from_fixture_structure(self):
        page, errors = self._open()
        page.fill("#editor", f"actor p = protein from={TEST_ID} chain=A\n0-2s p rotate")
        self.assertTrue(_poll(page, "() => window.__studio.ready()"), "scene never became ready")
        diag = page.locator("#diagnostics").inner_text()
        self.assertNotIn("generic shape", diag)
        self.assertEqual(page.evaluate("() => window.__studio.structureCount()"), 1)
        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()

    def test_export_downloads_a_webm(self):
        page, errors = self._open()
        options = page.eval_on_selector_all("#export option", "os => os.map(o => o.value)")
        if "webm:720p" not in options:
            self.skipTest("this Chromium build cannot encode WebM")
        self._load_example(page)
        with page.expect_download(timeout=300000) as dl_info:
            page.select_option("#export", "webm:720p")
        path = dl_info.value.path()
        with open(path, "rb") as fh:
            data = fh.read()
        self.assertTrue(dl_info.value.suggested_filename.endswith("-720p.webm"))
        self.assertEqual(data[:4], bytes.fromhex("1A45DFA3"), "not an EBML/WebM file")
        self.assertGreater(len(data), 50_000)
        expect(page.locator("#exportStatus")).to_contain_text("Saved", timeout=5000)
        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()

    def test_editor_starts_empty(self):
        page, errors = self._open()
        self.assertEqual(page.locator("#editor").input_value(), "")
        self.assertEqual(errors, [], f"uncaught JS errors: {errors}")
        page.close()


if __name__ == "__main__":
    unittest.main()
