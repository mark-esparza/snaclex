"""Tests for LIGSITE pocket detection and its scoring helpers."""

import unittest

from snaclex import pockets
from tests.fixtures import atom, structure


def _lining(res_names):
    return [
        {"res_id": f"A/{rn}{i}", "res_name": rn, "res_seq": i, "chain": "A"}
        for i, rn in enumerate(res_names)
    ]


class TestPocketHelpers(unittest.TestCase):
    def test_detect_pockets_too_few_atoms(self):
        s = structure([atom("C", i, 0, 0, res_seq=i) for i in range(5)])
        self.assertEqual(pockets.detect_pockets(s), [])

    def test_detect_pockets_returns_list(self):
        # A compact blob: should run end-to-end and return a (possibly empty) list.
        atoms = []
        n = 0
        for x in range(-2, 3):
            for y in range(-2, 3):
                atoms.append(atom("C", float(x), float(y), 0.0, res_seq=n))
                n += 1
        res = pockets.detect_pockets(structure(atoms))
        self.assertIsInstance(res, list)

    def test_assess_druggable(self):
        a = pockets._assess(volume=300, mean_psp=7.0, lining=_lining(["LEU"] * 10))
        self.assertEqual(a["tier"], "druggable")
        self.assertTrue(0 <= a["druggability_score"] <= 100)
        for key in ("volume", "enclosure", "hydrophobicity", "polarity", "aromaticity"):
            self.assertIn(key, a["subscores"])

    def test_assess_bare_pocket(self):
        a = pockets._assess(volume=10, mean_psp=5.0, lining=[])
        self.assertEqual(a["tier"], "pocket")
        self.assertTrue(0 <= a["druggability_score"] <= 100)

    def test_assess_score_bounds(self):
        # Even pathological inputs must clamp into [0, 100].
        a = pockets._assess(volume=99999, mean_psp=7.0, lining=_lining(["LEU"] * 50))
        self.assertLessEqual(a["druggability_score"], 100)

    def test_lining_residues(self):
        cell = pockets._CellGrid(
            [atom("C", 1.0, 0, 0, name="CB", res_name="LEU", res_seq=5, chain="A")],
            cell=5.0,
        )
        out = pockets._lining_residues([(0.0, 0.0, 0.0)], cell, cutoff=4.0)
        self.assertEqual(out[0]["res_id"], "A/LEU5")
        self.assertNotIn("_d2", out[0])  # internal field must be stripped


if __name__ == "__main__":
    unittest.main()
