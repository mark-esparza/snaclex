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


NUCLEIC_PDB = """\
ATOM      1  N   ALA A  10      11.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A  10      12.000  10.000  10.000  1.00 20.00           C
ATOM      3  P    DA B   1       5.000   5.000   5.000  1.00 30.00           P
ATOM      4  OP1  DA B   1       6.000   5.000   5.000  1.00 30.00           O
ATOM      5  N1   DA B   1       7.000   5.000   5.000  1.00 30.00           N
ATOM      6  P    DT B   2       8.000   6.000   5.000  1.00 30.00           P
END
"""


class TestNucleicParsing(unittest.TestCase):
    def setUp(self):
        self.s = pdbparse.parse_pdb(NUCLEIC_PDB)

    def test_nucleotides_classified(self):
        # Two protein atoms (ALA), four nucleic atoms (DA + DT).
        self.assertEqual(len(self.s.protein_atoms), 2)
        self.assertEqual(len(self.s.nucleic_atoms), 4)

    def test_nucleic_chains(self):
        self.assertEqual(self.s.nucleic_chains, ["B"])
        # Nucleotide chain B is not a protein chain.
        self.assertEqual(self.s.chains, ["A"])

    def test_nucleotides_not_components(self):
        # ATOM-record nucleotides must not be treated as hetero ligands.
        self.assertEqual(self.s.components, [])

    def test_default_structure_has_empty_nucleic(self):
        # Backward-compatible defaults for code that builds Structure directly.
        s = pdbparse.Structure(atoms=[], protein_atoms=[], components=[], chains=[])
        self.assertEqual(s.nucleic_atoms, [])
        self.assertEqual(s.nucleic_chains, [])


class TestParseMmcif(unittest.TestCase):
    def setUp(self):
        from tests.fixtures import mmcif_text
        rows = [
            {"element": "N", "name": "N", "comp": "MET", "seq": 1, "x": 1.0},
            {"element": "C", "name": "CA", "comp": "MET", "seq": 1, "x": 2.0},
            # An alternate location for the same atom — only the first is kept.
            {"element": "C", "name": "CB", "comp": "MET", "seq": 1, "x": 3.0, "alt": "A"},
            {"element": "C", "name": "CB", "comp": "MET", "seq": 1, "x": 3.5, "alt": "B"},
            {"element": "O", "name": "O", "comp": "MET", "seq": 1, "x": 4.0},
            # A hetero ligand component.
            {"group": "HETATM", "element": "C", "name": "C1", "comp": "LIG",
             "chain": "B", "seq": 900, "x": 20.0},
            {"group": "HETATM", "element": "O", "name": "O1", "comp": "LIG",
             "chain": "B", "seq": 900, "x": 21.0},
            # A second model — must be ignored.
            {"element": "N", "name": "N", "comp": "GLY", "chain": "A", "seq": 2,
             "x": 99.0, "model": 2},
        ]
        self.text = mmcif_text(rows)

    def test_auto_detected_as_mmcif(self):
        s = pdbparse.parse_structure(self.text)
        # 4 protein atoms (one CB altloc dropped), GLY from model 2 excluded.
        self.assertEqual(len(s.protein_atoms), 4)

    def test_components_and_chains(self):
        s = pdbparse.parse_mmcif(self.text)
        names = {c.res_name for c in s.components}
        self.assertIn("LIG", names)
        self.assertEqual(s.chains, ["A"])
        lig = next(c for c in s.components if c.res_name == "LIG")
        self.assertEqual(lig.kind, "ligand")

    def test_coordinates_parsed(self):
        s = pdbparse.parse_mmcif(self.text)
        n_atom = next(a for a in s.protein_atoms if a.name == "N")
        self.assertAlmostEqual(n_atom.x, 1.0)

    def test_subset_chain(self):
        # Build a 2-chain structure and confirm subsetting keeps only one chain.
        from tests.fixtures import mmcif_text
        rows = (
            [{"element": "C", "name": f"A{i}", "chain": "A", "seq": i} for i in range(4)]
            + [{"element": "C", "name": f"B{i}", "chain": "B", "seq": i} for i in range(3)]
        )
        s = pdbparse.parse_mmcif(mmcif_text(rows))
        sub = pdbparse.subset_chain(s, "A")
        self.assertEqual(sub.chains, ["A"])
        self.assertEqual(len(sub.protein_atoms), 4)
        self.assertTrue(all(a.chain == "A" for a in sub.atoms))

    def test_to_pdb_roundtrip(self):
        # mmCIF -> Structure -> PDB text -> Structure preserves atoms/coords, so
        # the 3Dmol viewer (which reads "pdb") can render an mmCIF-sourced entry.
        s1 = pdbparse.parse_mmcif(self.text)
        pdb_text = pdbparse.to_pdb(s1)
        s2 = pdbparse.parse_pdb(pdb_text)
        self.assertEqual(len(s2.protein_atoms), len(s1.protein_atoms))
        self.assertEqual(len(s2.atoms), len(s1.atoms))
        n_atom = next(a for a in s2.protein_atoms if a.name == "N")
        self.assertAlmostEqual(n_atom.x, 1.0)
        self.assertTrue(pdb_text.strip().endswith("END"))


if __name__ == "__main__":
    unittest.main()
