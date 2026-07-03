"""Offline tests for the ANM normal-mode motion module.

The eigensolver is cross-checked against a brute-force full Jacobi solve, and the
end-to-end analysis is checked on a compact (space-filling) synthetic fold.
"""

import math
import random
import unittest

from snaclex import motion
from snaclex.pdbparse import Atom
from tests.fixtures import structure as make_structure


def _compact_blob(n=180, seed=5, two_domain=True):
    random.seed(seed)
    atoms = []
    i = 0
    domains = [(0, n // 2), (22, n - n // 2)] if two_domain else [(0, n)]
    for cx, count in domains:
        placed = 0
        while placed < count:
            x = random.uniform(-11, 11)
            y = random.uniform(-11, 11)
            z = random.uniform(-11, 11)
            if x * x + y * y + z * z <= 121:
                atoms.append(Atom(i + 1, "CA", "ALA", "A", i + 1, "",
                                  cx + x, y, z, "C", False))
                i += 1
                placed += 1
    return make_structure(protein=atoms)


class TestEigensolver(unittest.TestCase):
    def test_subspace_matches_full_solve(self):
        # Small compact bead set: solver's non-trivial eigenvalues must match the
        # brute-force full symmetric eigensolve.
        random.seed(2)
        pts = [(random.uniform(0, 10), random.uniform(0, 10), random.uniform(0, 10))
               for _ in range(12)]
        H, dim = motion._hessian(pts, cutoff=8.0, gamma=1.0)
        full = sorted(motion._jacobi(H, dim, sweeps=100)[0])[6:12]
        low = [lam for lam, _ in motion._lowest_modes(H, dim, pts, 6)]
        for a, b in zip(full, low):
            self.assertAlmostEqual(a, b, places=3)

    def test_rigid_modes_are_null_space(self):
        # Rigid-body modes must satisfy H·r ≈ 0.
        random.seed(4)
        pts = [(random.uniform(0, 12), random.uniform(0, 12), random.uniform(0, 12))
               for _ in range(15)]
        H, dim = motion._hessian(pts, cutoff=9.0, gamma=1.0)
        for r in motion._rigid_basis(pts):
            Hr = motion._matvec(H, dim, r)
            self.assertLess(math.sqrt(sum(x * x for x in Hr)), 1e-6)


class TestAnalyze(unittest.TestCase):
    def test_modes_ascending_and_separated(self):
        res = motion.analyze(_compact_blob())
        self.assertTrue(res["available"])
        self.assertEqual(res["n_modes"], motion.N_MODES_DEFAULT)
        eigs = [m["eigenvalue"] for m in res["modes"]]
        # Ascending, non-negative, and the first internal mode is clearly > 0.
        self.assertTrue(all(eigs[i] <= eigs[i + 1] + 1e-9 for i in range(len(eigs) - 1)))
        self.assertTrue(all(e >= 0 for e in eigs))
        self.assertGreater(eigs[0], 1e-3)

    def test_bead_disp_normalized_and_shaped(self):
        res = motion.analyze(_compact_blob())
        self.assertEqual(len(res["bead_members"]), res["n_beads"])
        for m in res["modes"]:
            self.assertEqual(len(m["bead_disp"]), res["n_beads"])
            maxmag = max(math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2) for v in m["bead_disp"])
            self.assertAlmostEqual(maxmag, 1.0, places=3)
            self.assertTrue(0.0 <= m["collectivity"] <= 1.0)

    def test_coarse_graining_cap(self):
        res = motion.analyze(_compact_blob(n=300))
        self.assertLessEqual(res["n_beads"], motion.MAX_BEADS)
        self.assertTrue(res["coarse_grained"])
        # Every residue is represented by exactly one bead.
        total = sum(len(g) for g in res["bead_members"])
        self.assertEqual(total, res["residue_count"])

    def test_too_small_unavailable(self):
        atoms = [Atom(i + 1, "CA", "ALA", "A", i + 1, "", i * 3.8, 0, 0, "C", False)
                 for i in range(4)]
        res = motion.analyze(make_structure(protein=atoms))
        self.assertFalse(res["available"])


if __name__ == "__main__":
    unittest.main()
