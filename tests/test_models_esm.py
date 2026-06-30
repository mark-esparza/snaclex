"""Tests for the ESM model integration (offline; injectable transport seam)."""

import os
import unittest
from unittest import mock

from snaclex import models_esm


class TestAvailability(unittest.TestCase):
    def test_unavailable_without_token(self):
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertFalse(models_esm.available())
            # Public calls degrade gracefully (no raise) when no token/seam.
            self.assertFalse(models_esm.score_variants("ACDEF", ["A1V"])["available"])
            self.assertFalse(models_esm.fold_sequence("ACDEF")["available"])

    def test_available_with_token(self):
        with mock.patch.dict(os.environ, {"ESM_API_KEY": "x"}, clear=True):
            self.assertTrue(models_esm.available())
            self.assertTrue(models_esm.config()["available"])


class TestScoreVariants(unittest.TestCase):
    def test_scoring_with_injected_scorer(self):
        seq = "ACDEFGHIKL"

        def fake(s, pos, wt, mut):
            return -7.0 if mut == "P" else 0.5

        out = models_esm.score_variants(seq, ["C2P", "D3E"], scorer=fake)
        self.assertTrue(out["available"])
        by_input = {r["input"]: r for r in out["scored"]}
        self.assertTrue(by_input["C2P"]["scored"])
        self.assertEqual(by_input["C2P"]["effect"], "strongly disfavored")
        self.assertEqual(by_input["D3E"]["effect"], "tolerated")

    def test_wt_mismatch(self):
        out = models_esm.score_variants("ACDEF", ["A2V"], scorer=lambda *a: 0.0)
        r = out["scored"][0]
        self.assertFalse(r["scored"])
        self.assertIn("WT mismatch", r["reason"])

    def test_out_of_range(self):
        out = models_esm.score_variants("ACDEF", ["A99V"], scorer=lambda *a: 0.0)
        self.assertFalse(out["scored"][0]["scored"])

    def test_non_missense_skipped(self):
        out = models_esm.score_variants("ACDEF", ["A2*"], scorer=lambda *a: 0.0)
        self.assertFalse(out["scored"][0]["scored"])

    def test_scorer_failure_is_graceful(self):
        def boom(*a):
            raise RuntimeError("network down")

        out = models_esm.score_variants("ACDEF", ["A1V"], scorer=boom)
        self.assertFalse(out["scored"][0]["scored"])
        self.assertIn("unavailable", out["scored"][0]["reason"])


class TestFold(unittest.TestCase):
    def test_fold_with_injected_folder(self):
        out = models_esm.fold_sequence("ACDEFGHIKL", folder=lambda s: "ATOM ... \nEND")
        self.assertTrue(out["available"])
        self.assertEqual(out["source"], "esm3-predicted")
        self.assertIn("ATOM", out["pdb"])

    def test_fold_too_long(self):
        seq = "A" * (models_esm.MAX_FOLD_RESIDUES + 1)
        out = models_esm.fold_sequence(seq, folder=lambda s: "x")
        self.assertFalse(out["available"])
        self.assertIn("too long", out["reason"])

    def test_fold_failure_is_graceful(self):
        def boom(s):
            raise RuntimeError("no endpoint")

        out = models_esm.fold_sequence("ACDEF", folder=boom)
        self.assertFalse(out["available"])
        self.assertIn("unavailable", out["reason"])


if __name__ == "__main__":
    unittest.main()
