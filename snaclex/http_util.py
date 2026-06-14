"""Tiny HTTP helper around urllib with retry/backoff and friendly errors.

Public services like PubChem (PUG-REST) and RCSB rate-limit bursts and
occasionally return transient 5xx errors. Rather than failing the whole
analysis on a hiccup, transient failures are retried with exponential backoff;
genuine client errors (404 not-found, 400 bad-request) fail fast.
"""

from __future__ import annotations

import json
import os
import threading
import time
import urllib.error
import urllib.request

from . import cache as _cache

DEFAULT_TIMEOUT = 30
MAX_ATTEMPTS = 4
_BACKOFF = [0.6, 1.5, 3.0]  # seconds between attempts
_RETRY_CODES = {408, 429, 500, 502, 503, 504}
_UA = "SnaCleX/0.1 (research tool; +local)"

# Optional disk cache for GET responses, enabled by pointing SNACLEX_HTTP_CACHE
# at a writable directory (off by default, so dev/test behavior is unchanged).
# SNACLEX_HTTP_CACHE_TTL overrides the TTL in seconds.
_CACHE_LOCK = threading.Lock()
_CACHE = None
_CACHE_READY = False


def _get_cache():
    global _CACHE, _CACHE_READY
    if _CACHE_READY:
        return _CACHE
    with _CACHE_LOCK:
        if not _CACHE_READY:
            directory = os.environ.get("SNACLEX_HTTP_CACHE")
            if directory and directory.lower() not in ("0", "off", "false"):
                ttl = int(os.environ.get("SNACLEX_HTTP_CACHE_TTL") or 86400)
                try:
                    _CACHE = _cache.DiskCache(directory, ttl_seconds=ttl)
                except OSError:
                    _CACHE = None
            _CACHE_READY = True
    return _CACHE


def reset_cache():
    """Drop the cached singleton (used by tests that toggle the env var)."""
    global _CACHE, _CACHE_READY
    with _CACHE_LOCK:
        _CACHE = None
        _CACHE_READY = False


class FetchError(Exception):
    """Raised when an upstream fetch ultimately fails."""


class RateLimitError(FetchError):
    """Upstream is throttling us (HTTP 429/503); retried but still failing."""


def _read(url: str, timeout: int) -> bytes:
    """Fetch a URL with retry/backoff on transient errors (disk-cached if on)."""
    cache = _get_cache()
    if cache is not None:
        hit = cache.get(url)
        if hit is not None:
            return hit

    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    last: Exception | None = None
    throttled = False

    for attempt in range(MAX_ATTEMPTS):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            if cache is not None:
                cache.set(url, data)
            return data
        except urllib.error.HTTPError as exc:
            last = exc
            if exc.code in (429, 503):
                throttled = True
            if exc.code not in _RETRY_CODES:
                # Genuine client error (e.g. 404) — don't retry.
                raise FetchError(f"HTTP {exc.code} for {url}") from exc
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc

        if attempt < MAX_ATTEMPTS - 1:
            time.sleep(_BACKOFF[min(attempt, len(_BACKOFF) - 1)])

    if throttled:
        raise RateLimitError(
            "Upstream service is rate-limiting requests; please wait a few "
            "seconds and try again."
        ) from last
    reason = getattr(last, "reason", last)
    raise FetchError(f"Could not reach {url} after {MAX_ATTEMPTS} tries ({reason})") from last


def fetch_text(url: str, timeout: int = DEFAULT_TIMEOUT) -> str:
    return _read(url, timeout).decode("utf-8", errors="replace")


def fetch_bytes(url: str, timeout: int = DEFAULT_TIMEOUT) -> bytes:
    return _read(url, timeout)


def fetch_json(url: str, timeout: int = DEFAULT_TIMEOUT):
    raw = fetch_text(url, timeout=timeout)
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise FetchError(f"Invalid JSON from {url}") from exc
