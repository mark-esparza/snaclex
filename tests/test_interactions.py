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


class TestProteinInterface(unittest.TestCase):
    """profile_interface: contacts between two protein chain groups."""

    def test_salt_bridge_across_interface(self):
        # ARG NH1 on chain A near ASP OD1 on chain B -> salt bridge.
        a = [atom("N", 0, 0, 0, name="NH1", res_name="ARG", chain="A", res_seq=10)]
        b = [atom("O", 3.0, 0, 0, name="OD1", res_name="ASP", chain="B", res_seq=20)]
        prof = interactions.profile_interface(structure(a + b), ["A"], ["B"])
        self.assertEqual(prof["mode"], "protein-protein")
        self.assertGreaterEqual(prof["counts"]["salt_bridge"], 1)
        # Each side lists exactly its own residue.
        self.assertEqual(len(prof["interface_residues_a"]), 1)
        self.assertEqual(len(prof["interface_residues_b"]), 1)
        self.assertEqual(prof["interface_residues_a"][0]["chain"], "A")
        self.assertEqual(prof["interface_residues_b"][0]["chain"], "B")

    def test_hydrogen_bond_across_interface(self):
        a = [atom("N", 0, 0, 0, name="N", res_name="GLY", chain="A", res_seq=5)]
        b = [atom("O", 3.0, 0, 0, name="O", res_name="GLY", chain="B", res_seq=5)]
        prof = interactions.profile_interface(structure(a + b), ["A"], ["B"])
        self.assertGreaterEqual(prof["counts"]["hydrogen_bond"], 1)

    def test_no_self_contacts(self):
        # Two close residues on the SAME chain must not register as an interface.
        a = [
            atom("N", 0, 0, 0, name="NH1", res_name="ARG", chain="A", res_seq=10),
            atom("O", 3.0, 0, 0, name="OD1", res_name="ASP", chain="A", res_seq=11),
        ]
        prof = interactions.profile_interface(structure(a), ["A"], ["B"])
        self.assertEqual(prof["interaction_total"], 0)

    def test_aromatic_stacking_across_interface(self):
        # A PHE ring on chain A stacked ~4 A above a PHE ring on chain B.
        def phe_ring(chain, z):
            return [
                atom("C", 1.4, 0, z, name="CG", res_name="PHE", chain=chain, res_seq=1),
                atom("C", 0.7, 1.2, z, name="CD1", res_name="PHE", chain=chain, res_seq=1),
                atom("C", 0.7, -1.2, z, name="CD2", res_name="PHE", chain=chain, res_seq=1),
                atom("C", -0.7, 1.2, z, name="CE1", res_name="PHE", chain=chain, res_seq=1),
                atom("C", -0.7, -1.2, z, name="CE2", res_name="PHE", chain=chain, res_seq=1),
                atom("C", -1.4, 0, z, name="CZ", res_name="PHE", chain=chain, res_seq=1),
            ]
        atoms = phe_ring("A", 0.0) + phe_ring("B", 4.0)
        prof = interactions.profile_interface(structure(atoms), ["A"], ["B"])
        self.assertGreaterEqual(prof["counts"]["aromatic"], 1)

    def test_empty_when_group_missing(self):
        a = [atom("N", 0, 0, 0, name="N", res_name="GLY", chain="A", res_seq=5)]
        prof = interactions.profile_interface(structure(a), ["A"], ["Z"])
        self.assertEqual(prof["interaction_total"], 0)


class TestNucleicInterface(unittest.TestCase):
    """profile_nucleic_interface: protein .. DNA/RNA contacts."""

    def _nuc(self, element, x, y, z, *, name, res_name, res_seq):
        return atom(element, x, y, z, name=name, res_name=res_name,
                    chain="B", res_seq=res_seq)

    def test_phosphate_salt_bridge(self):
        prot = [atom("N", 0, 0, 0, name="NZ", res_name="LYS", chain="A", res_seq=10)]
        nuc = [self._nuc("O", 3.0, 0, 0, name="OP1", res_name="DA", res_seq=1)]
        prof = interactions.profile_nucleic_interface(structure(prot, nucleic=nuc))
        self.assertEqual(prof["mode"], "protein-nucleic")
        self.assertGreaterEqual(prof["counts"]["salt_bridge"], 1)
        self.assertEqual(prof["nucleic_chains"], ["B"])

    def test_hydrogen_bond_to_base(self):
        prot = [atom("O", 0, 0, 0, name="OG", res_name="SER", chain="A", res_seq=5)]
        nuc = [self._nuc("N", 3.0, 0, 0, name="N1", res_name="DA", res_seq=1)]
        prof = interactions.profile_nucleic_interface(structure(prot, nucleic=nuc))
        self.assertGreaterEqual(prof["counts"]["hydrogen_bond"], 1)
        self.assertEqual(len(prof["nucleic_residues"]), 1)

    def test_empty_without_nucleic(self):
        prot = [atom("N", 0, 0, 0, name="N", res_name="GLY", chain="A", res_seq=5)]
        prof = interactions.profile_nucleic_interface(structure(prot))
        self.assertEqual(prof["interaction_total"], 0)


if __name__ == "__main__":
    unittest.main()
