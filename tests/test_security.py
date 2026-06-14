"""Tests for the Phase 1 security/abuse controls in server.py.

These import ``server`` directly (importing does not start the HTTP server) and
exercise the rate limiter and input-sanitizer in isolation, plus the static
list of expensive endpoints and the CSP contents.
"""

import unittest

import server


class FakeClock:
    """A manually advanceable monotonic clock for deterministic rate tests."""

    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t

    def advance(self, dt):
        self.t += dt


class TestRateLimiter(unittest.TestCase):
    def test_allows_up_to_capacity_then_blocks(self):
        clock = FakeClock()
        rl = server.RateLimiter(rate=1.0, capacity=3, time_fn=clock)
        self.assertTrue(rl.allow("ip")[0])
        self.assertTrue(rl.allow("ip")[0])
        self.assertTrue(rl.allow("ip")[0])
        ok, retry = rl.allow("ip")
        self.assertFalse(ok)
        self.assertGreater(retry, 0)

    def test_refills_over_time(self):
        clock = FakeClock()
        rl = server.RateLimiter(rate=2.0, capacity=2, time_fn=clock)
        self.assertTrue(rl.allow("ip")[0])
        self.assertTrue(rl.allow("ip")[0])
        self.assertFalse(rl.allow("ip")[0])
        clock.advance(1.0)  # 2 tokens/sec -> 2 tokens back
        self.assertTrue(rl.allow("ip")[0])

    def test_keys_are_independent(self):
        clock = FakeClock()
        rl = server.RateLimiter(rate=1.0, capacity=1, time_fn=clock)
        self.assertTrue(rl.allow("a")[0])
        self.assertFalse(rl.allow("a")[0])
        self.assertTrue(rl.allow("b")[0])  # different IP, own bucket

    def test_eviction_bounds_memory(self):
        clock = FakeClock()
        rl = server.RateLimiter(rate=1.0, capacity=1, time_fn=clock, max_keys=8)
        for i in range(40):
            clock.advance(0.001)
            rl.allow(f"ip{i}")
        self.assertLessEqual(len(rl._buckets), 8)


class TestCleanText(unittest.TestCase):
    def test_strips_nul_and_control_chars(self):
        self.assertEqual(server.clean_text("as\x00pir\x07in"), "aspirin")

    def test_keeps_newlines_for_batch_input(self):
        self.assertEqual(server.clean_text("a\nb", max_len=50), "a\nb")

    def test_caps_length(self):
        self.assertEqual(len(server.clean_text("x" * 5000)), server.MAX_QUERY_LEN)

    def test_trims_whitespace(self):
        self.assertEqual(server.clean_text("  aspirin  "), "aspirin")

    def test_handles_none(self):
        self.assertEqual(server.clean_text(None), "")


class TestSecurityConfig(unittest.TestCase):
    def test_expensive_endpoints_listed(self):
        # Heavy GET endpoints carry the stricter budget; dock/screen now run
        # through the async job queue (POST /api/jobs) instead.
        for ep in ("/api/pockets", "/api/evolution"):
            self.assertIn(ep, server.EXPENSIVE_ENDPOINTS)
        self.assertNotIn("/api/dock", server.EXPENSIVE_ENDPOINTS)

    def test_csp_allows_required_sources(self):
        self.assertIn("https://3Dmol.org", server._CSP)
        self.assertIn("https://pubchem.ncbi.nlm.nih.gov", server._CSP)
        self.assertIn("frame-ancestors 'none'", server._CSP)
        self.assertIn("object-src 'none'", server._CSP)


if __name__ == "__main__":
    unittest.main()
