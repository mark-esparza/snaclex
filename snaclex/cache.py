"""TTL'd, size-bounded disk cache for upstream responses (pure stdlib).

Public scientific APIs (RCSB, PubChem, ChEMBL, Pfam) are slow and rate-limited,
and an in-memory cache is lost whenever the process restarts (frequent on
free-tier hosts). This stores raw response bytes on disk keyed by URL, with an
expiry stamp, so repeat lookups are instant and upstream load is reduced across
restarts.

Only opaque bytes are stored (no pickling of Python objects), so reading the
cache can never execute code. Writes are atomic (temp file + ``os.replace``).
"""

from __future__ import annotations

import hashlib
import os
import struct
import threading
import time

_HEADER = struct.Struct(">d")  # 8-byte big-endian expiry timestamp


class DiskCache:
    def __init__(self, directory, ttl_seconds=86400, max_entries=2000,
                 time_fn=time.time):
        self.dir = directory
        self.ttl = ttl_seconds
        self.max_entries = max_entries
        self._time = time_fn
        self._lock = threading.Lock()
        os.makedirs(self.dir, exist_ok=True)

    def _path(self, key: str) -> str:
        h = hashlib.sha256(key.encode("utf-8")).hexdigest()
        return os.path.join(self.dir, h + ".bin")

    def get(self, key: str):
        path = self._path(key)
        try:
            with open(path, "rb") as fh:
                head = fh.read(_HEADER.size)
                if len(head) < _HEADER.size:
                    return None
                (expiry,) = _HEADER.unpack(head)
                if self._time() > expiry:
                    self._safe_remove(path)
                    return None
                return fh.read()
        except (FileNotFoundError, OSError):
            return None

    def set(self, key: str, value: bytes):
        path = self._path(key)
        expiry = self._time() + self.ttl
        tmp = f"{path}.{os.getpid()}.{threading.get_ident()}.tmp"
        try:
            with open(tmp, "wb") as fh:
                fh.write(_HEADER.pack(expiry))
                fh.write(value)
            os.replace(tmp, path)
        except OSError:
            self._safe_remove(tmp)
            return
        self._evict_if_needed()

    def _evict_if_needed(self):
        with self._lock:
            try:
                entries = [e for e in os.scandir(self.dir) if e.name.endswith(".bin")]
            except OSError:
                return
            if len(entries) <= self.max_entries:
                return
            entries.sort(key=lambda e: e.stat().st_mtime)
            for e in entries[: len(entries) - self.max_entries]:
                self._safe_remove(e.path)

    @staticmethod
    def _safe_remove(path):
        try:
            os.remove(path)
        except OSError:
            pass
