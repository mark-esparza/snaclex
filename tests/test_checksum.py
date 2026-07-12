"""Tests for sequence checksums (the identity bridge). Offline."""

import unittest

from snaclex import checksum


class TestChecksums(unittest.TestCase):
    def test_deterministic(self):
        a = checksum.checksums("ACDEFGHIK")
        b = checksum.checksums("ACDEFGHIK")
        self.assertEqual(a, b)

    def test_different_sequences_differ(self):
        self.assertNotEqual(checksum.crc64("ACDEF"), checksum.crc64("ACDEG"))

    def test_crc64_is_16_hex(self):
        c = checksum.crc64("MEEPQSDPSV")
        self.assertEqual(len(c), 16)
        int(c, 16)  # parses as hex

    def test_whitespace_and_case_normalized(self):
        a = checksum.checksums("acd efg")
        b = checksum.checksums("ACDEFG")
        self.assertEqual(a["crc64"], b["crc64"])

    def test_md5_matches_hashlib(self):
        import hashlib
        self.assertEqual(checksum.md5("ACDEF"),
                         hashlib.md5(b"ACDEF").hexdigest())


if __name__ == "__main__":
    unittest.main()
