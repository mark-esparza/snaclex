"""Tests for the DiskCache and its wiring into http_util."""

import os
import tempfile
import unittest
from unittest import mock

from snaclex import http_util
from snaclex.cache import DiskCache


class FakeClock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t

    def advance(self, dt):
        self.t += dt


class TestDiskCache(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="snaclex-cache-test-")

    def test_set_get_roundtrip(self):
        c = DiskCache(self.dir)
        c.set("http://x/1", b"payload")
        self.assertEqual(c.get("http://x/1"), b"payload")

    def test_miss_returns_none(self):
        self.assertIsNone(DiskCache(self.dir).get("http://x/missing"))

    def test_ttl_expiry(self):
        clock = FakeClock()
        c = DiskCache(self.dir, ttl_seconds=10, time_fn=clock)
        c.set("k", b"v")
        self.assertEqual(c.get("k"), b"v")
        clock.advance(11)
        self.assertIsNone(c.get("k"))  # expired

    def test_eviction_bounds_entries(self):
        c = DiskCache(self.dir, max_entries=5)
        for i in range(20):
            c.set(f"k{i}", b"v")
        count = len([n for n in os.listdir(self.dir) if n.endswith(".bin")])
        self.assertLessEqual(count, 5)


class FakeResp:
    def __init__(self, body):
        self._b = body

    def read(self):
        return self._b

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


class TestHttpUtilCaching(unittest.TestCase):
    """Caching is off by default; when enabled, a second fetch is served locally."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="snaclex-http-cache-test-")
        http_util.reset_cache()
        self.addCleanup(http_util.reset_cache)

    def test_off_by_default(self):
        # No SNACLEX_HTTP_CACHE -> every call hits the network.
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("SNACLEX_HTTP_CACHE", None)
            http_util.reset_cache()
            with mock.patch.object(
                http_util.urllib.request, "urlopen",
                side_effect=[FakeResp(b"a"), FakeResp(b"b")],
            ) as m:
                self.assertEqual(http_util.fetch_text("http://x"), "a")
                self.assertEqual(http_util.fetch_text("http://x"), "b")
                self.assertEqual(m.call_count, 2)

    def test_second_fetch_served_from_cache(self):
        with mock.patch.dict(os.environ, {"SNACLEX_HTTP_CACHE": self.dir}):
            http_util.reset_cache()
            with mock.patch.object(
                http_util.urllib.request, "urlopen",
                side_effect=[FakeResp(b"once")],
            ) as m:
                self.assertEqual(http_util.fetch_text("http://x"), "once")
                self.assertEqual(http_util.fetch_text("http://x"), "once")
                self.assertEqual(m.call_count, 1)  # second call cached


if __name__ == "__main__":
    unittest.main()
