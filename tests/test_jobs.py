"""Tests for the in-process async JobManager."""

import time
import unittest

from snaclex.jobs import JobManager


def _wait(jm, job_id, timeout=5.0):
    """Block until a job leaves queued/running (real worker threads)."""
    end = time.time() + timeout
    while time.time() < end:
        st = jm.status(job_id)
        if st and st["status"] in ("done", "error"):
            return st
        time.sleep(0.01)
    raise AssertionError(f"job {job_id} did not finish in time")


class FakeClock:
    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t

    def advance(self, dt):
        self.t += dt


class TestJobManager(unittest.TestCase):
    def test_successful_job(self):
        jm = JobManager(max_workers=2)
        self.addCleanup(jm.shutdown)
        jid = jm.submit(lambda a, b: a + b, 2, 3)
        st = _wait(jm, jid)
        self.assertEqual(st["status"], "done")
        self.assertEqual(st["result"], 5)
        self.assertIsNone(st["error"])

    def test_failed_job_captures_message(self):
        jm = JobManager(max_workers=2)
        self.addCleanup(jm.shutdown)

        def boom():
            raise ValueError("kaboom")

        st = _wait(jm, jm.submit(boom))
        self.assertEqual(st["status"], "error")
        self.assertIn("kaboom", st["error"])

    def test_unknown_job_returns_none(self):
        jm = JobManager(max_workers=1)
        self.addCleanup(jm.shutdown)
        self.assertIsNone(jm.status("does-not-exist"))

    def test_status_is_a_copy(self):
        jm = JobManager(max_workers=1)
        self.addCleanup(jm.shutdown)
        jid = jm.submit(lambda: 1)
        _wait(jm, jid)
        snap = jm.status(jid)
        snap["status"] = "tampered"
        self.assertEqual(jm.status(jid)["status"], "done")

    def test_finished_jobs_are_gc_d_after_ttl(self):
        clock = FakeClock()
        jm = JobManager(max_workers=1, ttl_seconds=100, time_fn=clock)
        self.addCleanup(jm.shutdown)
        jid = jm.submit(lambda: 1)
        _wait(jm, jid)
        self.assertIsNotNone(jm.status(jid))
        clock.advance(200)              # age the finished job past its TTL
        jm.submit(lambda: 2)            # submitting triggers GC
        self.assertIsNone(jm.status(jid))


if __name__ == "__main__":
    unittest.main()
