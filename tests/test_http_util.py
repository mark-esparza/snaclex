"""Tests for the HTTP fetch layer's retry/backoff behavior.

This is the single network choke point for every upstream client (RCSB,
PubChem, ChEMBL, Pfam). By monkeypatching ``urlopen`` here we establish the seam
that lets the higher-level clients be tested fully offline.
"""

import unittest
import urllib.error
from unittest import mock

from snaclex import http_util
from snaclex.http_util import FetchError, RateLimitError


class FakeResp:
    def __init__(self, body: bytes):
        self._body = body

    def read(self):
        return self._body

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _http_error(code):
    return urllib.error.HTTPError("http://x", code, f"err{code}", None, None)


class TestFetch(unittest.TestCase):
    def setUp(self):
        # Never actually sleep during retry tests.
        patcher = mock.patch.object(http_util.time, "sleep", lambda *a: None)
        patcher.start()
        self.addCleanup(patcher.stop)

    def _patch_urlopen(self, side_effect):
        return mock.patch.object(
            http_util.urllib.request, "urlopen", side_effect=side_effect
        )

    def test_success(self):
        with self._patch_urlopen([FakeResp(b"hello")]) as m:
            self.assertEqual(http_util.fetch_text("http://x"), "hello")
            self.assertEqual(m.call_count, 1)

    def test_404_fails_fast_without_retry(self):
        with self._patch_urlopen([_http_error(404)]) as m:
            with self.assertRaises(FetchError):
                http_util.fetch_text("http://x")
            self.assertEqual(m.call_count, 1)  # genuine client error: no retry

    def test_retries_then_succeeds(self):
        seq = [_http_error(503), FakeResp(b"ok")]
        with self._patch_urlopen(seq) as m:
            self.assertEqual(http_util.fetch_text("http://x"), "ok")
            self.assertEqual(m.call_count, 2)

    def test_persistent_throttle_raises_ratelimit(self):
        seq = [_http_error(429)] * http_util.MAX_ATTEMPTS
        with self._patch_urlopen(seq) as m:
            with self.assertRaises(RateLimitError):
                http_util.fetch_text("http://x")
            self.assertEqual(m.call_count, http_util.MAX_ATTEMPTS)

    def test_fetch_json_parses(self):
        with self._patch_urlopen([FakeResp(b'{"a": 1}')]):
            self.assertEqual(http_util.fetch_json("http://x"), {"a": 1})

    def test_fetch_json_invalid_raises(self):
        with self._patch_urlopen([FakeResp(b"not json")]):
            with self.assertRaises(FetchError):
                http_util.fetch_json("http://x")


if __name__ == "__main__":
    unittest.main()
