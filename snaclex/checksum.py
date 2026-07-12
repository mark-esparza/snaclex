"""Sequence checksums — the identity key for the cross-reference graph.

Two proteins are the *same sequence* iff their checksums match; they are the
*same protein record* only when checksum **and** taxon **and** isoform role match
(see docs/platform/05-identifier-resolution.md). This lets a raw FASTA paste
resolve to a UniProt/UniParc record without a name lookup, and stops orthologs
(same-ish sequence, different taxon) from being merged.

* ``crc64`` implements the SWISS-PROT / UniParc CRC-64 (ISO 3309, POLY64REV
  0xd8000000), ported from the canonical BioPython implementation so it matches
  UniProt's ``sequence.crc64`` field. Used to query UniParc ``checksum:``.
* ``md5``/``sha256`` are exact stdlib digests, used as the primary internal key.
"""

from __future__ import annotations

import hashlib

_TABLE_H: list[int] = []
_TABLE_L: list[int] = []


def _init_tables() -> None:
    for i in range(256):
        low = i
        part_h = 0
        for _ in range(8):
            rflag = low & 1
            low >>= 1
            if part_h & 1:
                low |= 1 << 31
            part_h >>= 1
            if rflag:
                part_h ^= 0xD8000000
        _TABLE_H.append(part_h)
        _TABLE_L.append(low)


_init_tables()


def crc64(seq: str) -> str:
    """Return the SWISS-PROT/UniParc CRC-64 as a 16-char uppercase hex string.

    Matches UniProt's ``sequence.crc64`` value (without the ``CRC64`` prefix).
    """
    crch = 0
    crcl = 0
    for ch in seq:
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        idx = (crcl ^ ord(ch)) & 0xFF
        crch = temp1h ^ _TABLE_H[idx]
        crcl = temp1l ^ _TABLE_L[idx]
    return "%08X%08X" % (crch, crcl)


def md5(seq: str) -> str:
    return hashlib.md5(seq.encode("ascii", "replace")).hexdigest()


def sha256(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii", "replace")).hexdigest()


def checksums(seq: str) -> dict:
    """All three checksums for a sequence (uppercased, whitespace-stripped)."""
    norm = "".join(seq.split()).upper()
    return {"crc64": crc64(norm), "md5": md5(norm), "sha256": sha256(norm)}
