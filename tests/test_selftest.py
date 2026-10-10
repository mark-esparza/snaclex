"""Tests for the upstream connectivity self-test.

The checks themselves need the network, so what is tested here is the runner:
that a down upstream is reported rather than raised, that the report shape is
stable, and that the fail-fast fetch policy is thread-local (the server is
threaded, so a diagnostic must not alter concurrent real requests).
"""

import threading
import unittest
import urllib.error

from snaclex import http_util, selftest


def _ok(detail="fine"):
    return lambda: detail


def _boom(exc=RuntimeError("upstream down")):
    def fail():
        raise exc
    return fail


def _checks(*specs):
    """Build CHECKS-shaped tuples: (id, source, host, proves, callable)."""
    return [
        (f"c{i}", src, "example.org", "proves something", fn)
        for i, (src, fn) in enumerate(specs)
    ]


class TestRunChecks(unittest.TestCase):
    def test_all_passing(self):
        r = selftest.run_checks(_checks(("RCSB", _ok("a")), ("PubChem", _ok("b"))))
        self.assertTrue(r["all_ok"])
        self.assertEqual((r["n_passed"], r["n_failed"], r["n_checks"]), (2, 0, 2))
        self.assertEqual([c["detail"] for c in r["checks"]], ["a", "b"])

    def test_a_down_upstream_is_reported_not_raised(self):
        r = selftest.run_checks(_checks(("RCSB", _ok()), ("PubChem", _boom())))
        self.assertFalse(r["all_ok"])
        self.assertEqual(r["n_failed"], 1)
        bad = [c for c in r["checks"] if not c["ok"]][0]
        self.assertIn("upstream down", bad["error"])
        self.assertIn("RuntimeError", bad["error"])

    def test_one_failure_does_not_stop_later_checks(self):
        r = selftest.run_checks(
            _checks(("A", _boom()), ("B", _ok("ran")), ("C", _ok("ran")))
        )
        self.assertEqual(r["n_passed"], 2)
        self.assertEqual(len(r["checks"]), 3)

    def test_per_source_rollup(self):
        r = selftest.run_checks(
            _checks(("RCSB", _ok()), ("RCSB", _boom()), ("PubChem", _ok()))
        )
        self.assertEqual(r["by_source"]["RCSB"], {"ok": 1, "failed": 1})
        self.assertEqual(r["by_source"]["PubChem"], {"ok": 1, "failed": 0})

    def test_latency_recorded(self):
        ticks = iter([0.0, 0.25, 0.25, 0.5])
        r = selftest.run_checks(
            _checks(("A", _ok()), ("B", _ok())), timer=lambda: next(ticks)
        )
        self.assertEqual([c["ms"] for c in r["checks"]], [250, 250])

    def test_report_is_json_serializable(self):
        import json
        json.dumps(selftest.run_checks(_checks(("A", _ok()), ("B", _boom()))))

    def test_real_check_table_covers_every_integrated_source(self):
        sources = {c[1] for c in selftest.CHECKS}
        self.assertEqual(
            sources, {"RCSB PDB", "PubChem", "ChEMBL", "InterPro / Pfam"}
        )
        hosts = {c[2] for c in selftest.CHECKS}
        for host in ("files.rcsb.org", "data.rcsb.org", "search.rcsb.org",
                     "pubchem.ncbi.nlm.nih.gov", "www.ebi.ac.uk"):
            self.assertIn(host, hosts)

    def test_formatted_report_names_failures(self):
        r = selftest.run_checks(_checks(("RCSB", _boom())))
        text = selftest.format_report(r)
        self.assertIn("FAIL", text)
        self.assertIn("0/1 checks passed", text)
        self.assertIn("egress", text)


class TestFailFastPolicy(unittest.TestCase):
    def test_defaults_are_unchanged_outside_the_context(self):
        self.assertIsNone(getattr(http_util._local, "attempts", None))

    def test_context_sets_and_restores(self):
        with http_util.fail_fast(attempts=2, timeout=5):
            self.assertEqual(http_util._local.attempts, 2)
            self.assertEqual(http_util._local.timeout, 5)
        self.assertIsNone(http_util._local.attempts)
        self.assertIsNone(http_util._local.timeout)

    def test_nested_contexts_restore_the_outer_value(self):
        with http_util.fail_fast(attempts=3):
            with http_util.fail_fast(attempts=1):
                self.assertEqual(http_util._local.attempts, 1)
            self.assertEqual(http_util._local.attempts, 3)

    def test_policy_does_not_leak_to_other_threads(self):
        """The server is threaded: a diagnostic must not change live requests."""
        seen = {}

        def worker():
            seen["attempts"] = getattr(http_util._local, "attempts", None)

        with http_util.fail_fast(attempts=1):
            t = threading.Thread(target=worker)
            t.start()
            t.join()
        self.assertIsNone(seen["attempts"])

    def test_retry_count_is_honoured(self):
        """One attempt means one call, and the message says so."""
        calls = []

        def fake_urlopen(req, timeout=None):
            calls.append(timeout)
            raise urllib.error.URLError("nope")

        import urllib.error
        import urllib.request
        real = urllib.request.urlopen
        urllib.request.urlopen = fake_urlopen
        try:
            with http_util.fail_fast(attempts=1, timeout=4):
                with self.assertRaises(http_util.FetchError) as ctx:
                    http_util.fetch_text("https://example.org/x")
        finally:
            urllib.request.urlopen = real
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0], 4)          # timeout clamped down
        self.assertIn("after 1 try", str(ctx.exception))

    def test_without_the_context_the_full_retry_budget_is_used(self):
        calls = []

        def fake_urlopen(req, timeout=None):
            calls.append(timeout)
            raise urllib.error.URLError("nope")

        import urllib.error
        import urllib.request
        real = urllib.request.urlopen
        urllib.request.urlopen = fake_urlopen
        real_backoff = http_util._BACKOFF
        http_util._BACKOFF = [0, 0, 0]          # keep the test fast
        try:
            with self.assertRaises(http_util.FetchError):
                http_util.fetch_text("https://example.org/x")
        finally:
            urllib.request.urlopen = real
            http_util._BACKOFF = real_backoff
        self.assertEqual(len(calls), http_util.MAX_ATTEMPTS)


if __name__ == "__main__":
    unittest.main()
