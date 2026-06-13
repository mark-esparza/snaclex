"""Tests for the pure-Python rigid-body docker."""

import unittest

from snaclex import docking
from snaclex.pdbparse import Component
from tests.fixtures import atom, structure


def _receptor():
    # A small cluster of protein atoms around the origin to give the grid signal.
    coords = [
        ("C", 3, 0, 0), ("O", -3, 0, 0), ("N", 0, 3, 0), ("C", 0, -3, 0),
        ("C", 2, 2, 1), ("O", -2, 2, -1), ("N", 2, -2, 1), ("C", -2, -2, -1),
        ("C", 4, 1, 0), ("C", -4, -1, 0), ("O", 1, 4, 0), ("N", -1, -4, 0),
    ]
    return structure(
        [atom(el, x, y, z, name=el, res_name="LEU", res_seq=i)
         for i, (el, x, y, z) in enumerate(coords)]
    )


class TestDocking(unittest.TestCase):
    def setUp(self):
        self.s = _receptor()
        self.lig = [
            {"element": "C", "x": 0.0, "y": 0.0, "z": 0.0},
            {"element": "O", "x": 1.4, "y": 0.0, "z": 0.0},
        ]

    def test_deterministic_with_seed(self):
        a = docking.dock(self.s, self.lig, (0, 0, 0), seeds=20, mc_steps=10, seed=5)
        b = docking.dock(self.s, self.lig, (0, 0, 0), seeds=20, mc_steps=10, seed=5)
        self.assertEqual(a["score"], b["score"])
        self.assertEqual(a["pose_coords"], b["pose_coords"])

    def test_output_schema(self):
        pose = docking.dock(self.s, self.lig, (0, 0, 0), seeds=20, mc_steps=10, seed=1)
        self.assertEqual(pose["n_heavy_atoms"], 2)
        self.assertEqual(len(pose["pose_coords"]), 2)
        self.assertEqual(pose["elements"], ["C", "O"])
        # LE is score-per-heavy-atom (computed from the unrounded score, so
        # compare against the 2-dp score with matching tolerance).
        self.assertAlmostEqual(pose["ligand_efficiency"], pose["score"] / 2, places=2)
        self.assertEqual(
            pose["search"], {"seeds": 20, "mc_steps": 10, "random_seed": 1}
        )

    def test_empty_ligand_raises(self):
        grid = docking.build_grid(self.s, (0, 0, 0))
        with self.assertRaises(ValueError):
            docking.dock_with_grid(grid, [], (0, 0, 0))

    def test_pose_to_component(self):
        pose = {"pose_coords": [(0, 0, 0), (1, 1, 1)], "elements": ["C", "O"]}
        comp = docking.pose_to_component(pose, "LIG")
        self.assertEqual(len(comp.atoms), 2)
        self.assertEqual(comp.chain, "X")
        self.assertTrue(all(a.is_hetero for a in comp.atoms))

    def test_pose_to_pdb(self):
        pose = {"pose_coords": [(0, 0, 0), (1, 1, 1)], "elements": ["C", "O"]}
        text = docking.pose_to_pdb(pose, "LIG")
        self.assertEqual(text.count("HETATM"), 2)
        self.assertTrue(text.strip().endswith("END"))

    def test_component_center(self):
        comp = Component("LIG", "X", 1, "", [
            atom("C", 0, 0, 0, hetero=True),
            atom("C", 2, 2, 2, hetero=True),
        ])
        self.assertEqual(docking.component_center(comp), (1.0, 1.0, 1.0))

    def test_rmsd_to_reference_perfect(self):
        ref = Component("LIG", "X", 1, "", [
            atom("C", 0, 0, 0, hetero=True),
            atom("O", 1, 1, 1, hetero=True),
        ])
        pose = {"pose_coords": [(0, 0, 0), (1, 1, 1)], "elements": ["C", "O"]}
        self.assertEqual(docking.rmsd_to_reference(pose, ref), 0.0)

    def test_rmsd_to_reference_no_element_match(self):
        ref = Component("LIG", "X", 1, "", [atom("C", 0, 0, 0, hetero=True)])
        pose = {"pose_coords": [(0, 0, 0)], "elements": ["N"]}
        self.assertIsNone(docking.rmsd_to_reference(pose, ref))


if __name__ == "__main__":
    unittest.main()
