"""Tests for the atomic interaction profiler.

Each test builds a minimal structure with geometry chosen to trigger exactly one
interaction class, so the geometric classifier can be checked in isolation.
"""

import unittest

from snaclex import interactions
from tests.fixtures import atom, component, structure


def _ligand(*atoms):
    return component("LIG", atoms, chain="X", res_seq=900)


class TestInteractions(unittest.TestCase):
    def test_hydrogen_bond(self):
        lig = _ligand(
            atom("O", 0, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
            atom("C", 1.5, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        )
        prot = [atom("N", 3.0, 0, 0, name="N", res_name="ALA", res_seq=10)]
        profile = interactions.profile_component(structure(prot, [lig]), lig)
        self.assertGreaterEqual(profile["counts"]["hydrogen_bond"], 1)

    def test_salt_bridge(self):
        # Ligand N near ASP carboxylate oxygen -> salt bridge (checked before H-bond).
        lig = _ligand(
            atom("N", 0, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
            atom("C", 1.5, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        )
        prot = [atom("O", 3.0, 0, 0, name="OD1", res_name="ASP", res_seq=20)]
        profile = interactions.profile_component(structure(prot, [lig]), lig)
        self.assertGreaterEqual(profile["counts"]["salt_bridge"], 1)
        self.assertEqual(profile["counts"]["hydrogen_bond"], 0)

    def test_hydrophobic(self):
        lig = _ligand(
            atom("C", 0, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
            atom("C", 1.5, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        )
        prot = [atom("C", 3.5, 0, 0, name="CB", res_name="LEU", res_seq=30)]
        profile = interactions.profile_component(structure(prot, [lig]), lig)
        self.assertGreaterEqual(profile["counts"]["hydrophobic"], 1)

    def test_metal_coordination(self):
        metal = component(
            "ZN",
            [atom("ZN", 0, 0, 0, name="ZN", res_name="ZN", chain="X",
                  res_seq=901, hetero=True)],
            chain="X",
            res_seq=901,
        )
        self.assertEqual(metal.kind, "metal")
        prot = [atom("O", 2.2, 0, 0, name="OD1", res_name="ASP", res_seq=40)]
        profile = interactions.profile_component(structure(prot, [metal]), metal)
        self.assertGreaterEqual(profile["counts"]["metal_coordination"], 1)

    def test_aromatic(self):
        # A PHE ring centered on the origin; a ligand atom sits above the centroid.
        ring = [
            atom("C", 1.4, 0, 0, name="CG", res_name="PHE", res_seq=50),
            atom("C", 0.7, 1.2, 0, name="CD1", res_name="PHE", res_seq=50),
            atom("C", 0.7, -1.2, 0, name="CD2", res_name="PHE", res_seq=50),
            atom("C", -0.7, 1.2, 0, name="CE1", res_name="PHE", res_seq=50),
            atom("C", -0.7, -1.2, 0, name="CE2", res_name="PHE", res_seq=50),
            atom("C", -1.4, 0, 0, name="CZ", res_name="PHE", res_seq=50),
        ]
        lig = _ligand(
            atom("C", 0, 0, 3.5, hetero=True, res_name="LIG", chain="X", res_seq=900),
            atom("N", 0, 0, 4.9, hetero=True, res_name="LIG", chain="X", res_seq=900),
        )
        profile = interactions.profile_component(structure(ring, [lig]), lig)
        self.assertGreaterEqual(profile["counts"]["aromatic"], 1)

    def test_profile_schema(self):
        lig = _ligand(
            atom("O", 0, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
            atom("C", 1.5, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        )
        prot = [atom("N", 3.0, 0, 0, name="N", res_name="ALA", res_seq=10)]
        profile = interactions.profile_component(structure(prot, [lig]), lig)
        for key in (
            "component", "interactions", "counts", "interaction_total",
            "contact_residues", "contact_residue_count",
        ):
            self.assertIn(key, profile)
        self.assertEqual(
            profile["interaction_total"], len(profile["interactions"])
        )
        # Per-residue summary types are JSON-serializable lists (not sets).
        for res in profile["contact_residues"]:
            self.assertIsInstance(res["types"], list)

    def test_no_contacts_when_far(self):
        lig = _ligand(
            atom("O", 0, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
            atom("C", 1.5, 0, 0, hetero=True, res_name="LIG", chain="X", res_seq=900),
        )
        prot = [atom("N", 50, 50, 50, name="N", res_name="ALA", res_seq=10)]
        profile = interactions.profile_component(structure(prot, [lig]), lig)
        self.assertEqual(profile["interaction_total"], 0)


if __name__ == "__main__":
    unittest.main()
