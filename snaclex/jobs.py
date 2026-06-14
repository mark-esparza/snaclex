"""In-process asynchronous job queue (pure stdlib).

Long-running analyses (docking, batch screening) are submitted here instead of
running inline in the HTTP request, so the connection isn't held open for the
whole compute (which trips proxy/browser timeouts on hosts like Render) and a
burst of work *queues* behind a bounded worker pool rather than being rejected.

Jobs are kept in memory with a TTL and garbage-collected once they finish and
age out. This is the dependency-free default; an optional Celery/Redis-backed
manager could be swapped in behind the same submit/status interface.
"""

from __future__ import annotations

import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor


class JobManager:
    """Submit callables to a bounded thread pool and poll their status.

    Status dicts have: ``status`` (queued|running|done|error), ``result``,
    ``error``, ``created`` and ``updated`` timestamps.
    """

    def __init__(self, max_workers=4, ttl_seconds=900, time_fn=time.monotonic):
        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        self._jobs: dict[str, dict] = {}
        self._lock = threading.Lock()
        self._ttl = ttl_seconds
        self._time = time_fn

    def submit(self, fn, *args, **kwargs) -> str:
        job_id = uuid.uuid4().hex
        now = self._time()
        with self._lock:
            self._jobs[job_id] = {
                "status": "queued",
                "result": None,
                "error": None,
                "created": now,
                "updated": now,
            }
        self._executor.submit(self._run, job_id, fn, args, kwargs)
        self._gc()
        return job_id

    def status(self, job_id: str):
        """Return a copy of the job's status dict, or None if unknown/expired."""
        with self._lock:
            job = self._jobs.get(job_id)
            return dict(job) if job is not None else None

    def _set(self, job_id, **fields):
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                return
            job.update(fields)
            job["updated"] = self._time()

    def _run(self, job_id, fn, args, kwargs):
        self._set(job_id, status="running")
        try:
            result = fn(*args, **kwargs)
            self._set(job_id, status="done", result=result)
        except Exception as exc:  # noqa: BLE001 - surfaced to the client as job error
            self._set(job_id, status="error", error=str(exc))

    def _gc(self):
        cutoff = self._time() - self._ttl
        with self._lock:
            stale = [
                k for k, v in self._jobs.items()
                if v["status"] in ("done", "error") and v["updated"] < cutoff
            ]
            for k in stale:
                self._jobs.pop(k, None)

    def shutdown(self):
        self._executor.shutdown(wait=False)
