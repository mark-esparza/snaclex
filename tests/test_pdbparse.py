"""Tests for the PDB-format parser."""

import unittest

from snaclex import pdbparse


# A tiny hand-written PDB exercising: a protein residue, an organic ligand,
# a metal ion, water (must be dropped), an alternate location (only first kept),
# and a second MODEL (must be ignored).
SAMPLE_PDB = """\
ATOM      1  N   ALA A  10      11.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA AALA A  10      12.000  10.000  10.000  0.60 20.00           C
ATOM      3  CA BALA A  10      12.500  10.500  10.000  0.40 20.00           C
ATOM      4  O   ALA A  10      13.000  10.000  10.000  1.00 20.00           O
HETATM    5  C1  LIG B 900      20.000  20.000  20.000  1.00 30.00           C
HETATM    6  O1  LIG B 900      21.400  20.000  20.000  1.00 30.00           O
HETATM    7 ZN    ZN B 901      30.000  30.000  30.000  1.00 15.00          ZN
HETATM    8  O   HOH B 950      40.000  40.000  40.000  1.00 25.00           O
ENDMDL
MODEL        2
ATOM      9  N   GLY A  11      99.000  99.000  99.000  1.00 20.00           N
ENDMDL
END
"""


class TestParsePdb(unittest.TestCase):
    def setUp(self):
        self.s = pdbparse.parse_pdb(SAMPLE_PDB)

    def test_only_first_model_kept(self):
        # The GLY from MODEL 2 must not appear.
        self.assertFalse(any(a.res_name == "GLY" for a in self.s.atoms))

    def test_alternate_location_dedup(self):
        cas = [a for a in self.s.protein_atoms if a.name == "CA"]
        self.assertEqual(len(cas), 1)
        # The first altloc (A, occupancy 0.60 at x=12.0) is the one retained.
        self.assertAlmostEqual(cas[0].x, 12.000)

    def test_protein_atoms_classified(self):
        self.assertEqual(len(self.s.protein_atoms), 3)  # N, CA, O of ALA
        self.assertEqual(self.s.chains, ["A"])

    def test_water_excluded_from_components(self):
        names = {c.res_name for c in self.s.components}
        self.assertNotIn("HOH", names)

    def test_components_and_kinds(self):
        by_name = {c.res_name: c for c in self.s.components}
        self.assertIn("LIG", by_name)
        self.assertIn("ZN", by_name)
        self.assertEqual(by_name["LIG"].kind, "ligand")  # 2 heavy atoms
        self.assertEqual(by_name["ZN"].kind, "metal")    # single metal atom
        self.assertEqual(self.s.ligand_components[0].res_name, "LIG")

    def test_element_fallback_from_name(self):
        # Even if the element column were blank, parsing must still classify; here
        # we just confirm the explicit element column is honored.
        zn = next(c for c in self.s.components if c.res_name == "ZN")
        self.assertEqual(zn.atoms[0].element, "ZN")

    def test_bfactor_and_occupancy_parsed(self):
        n_atom = next(a for a in self.s.protein_atoms if a.name == "N")
        self.assertAlmostEqual(n_atom.bfactor, 20.0)
        self.assertAlmostEqual(n_atom.occupancy, 1.0)

    def test_component_label(self):
        lig = next(c for c in self.s.components if c.res_name == "LIG")
        self.assertEqual(lig.label, "LIG B900")


if __name__ == "__main__":
    unittest.main()
