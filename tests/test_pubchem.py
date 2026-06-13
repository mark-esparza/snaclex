"""Tests for PubChem pure-Python helpers (no network).

The network-touching entry points (``lookup_compound``/``fetch_3d_atoms``) are
covered indirectly: their parsing/scoring helpers are unit-tested here, and the
shared fetch layer is tested in ``test_http_util``.
"""

import unittest

from snaclex import pubchem
from snaclex.http_util import FetchError


class TestResolveNamespace(unittest.TestCase):
    def test_numeric_is_cid(self):
        self.assertEqual(pubchem._resolve_namespace("2244"), ("cid", "2244"))

    def test_name_is_url_quoted(self):
        ns, value = pubchem._resolve_namespace("acetylsalicylic acid")
        self.assertEqual(ns, "name")
        self.assertEqual(value, "acetylsalicylic%20acid")

    def test_empty_raises(self):
        with self.assertRaises(FetchError):
            pubchem._resolve_namespace("  ")


class TestParseSdfAtoms(unittest.TestCase):
    SDF = (
        "\n"
        "  test\n"
        "\n"
        "  3  2  0  0  0  0            999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0\n"
        "    1.4000    0.0000    0.0000 O   0  0\n"
        "    2.0000    0.0000    0.0000 H   0  0\n"
        "M  END\n"
    )

    def test_parses_heavy_atoms_only(self):
        atoms = pubchem._parse_sdf_atoms(self.SDF)
        self.assertEqual(len(atoms), 2)  # H dropped
        self.assertEqual(atoms[0]["element"], "C")
        self.assertAlmostEqual(atoms[1]["x"], 1.4)

    def test_truncated_sdf_returns_empty(self):
        self.assertEqual(pubchem._parse_sdf_atoms("too\nshort\n"), [])


class TestToFloat(unittest.TestCase):
    def test_valid(self):
        self.assertEqual(pubchem._to_float("3.14"), 3.14)

    def test_invalid_returns_none(self):
        self.assertIsNone(pubchem._to_float("n/a"))
        self.assertIsNone(pubchem._to_float(None))


class TestLipinski(unittest.TestCase):
    def test_drug_like(self):
        props = {"molecular_weight": 180, "xlogp": 1.2,
                 "h_bond_donors": 1, "h_bond_acceptors": 4}
        out = pubchem._lipinski(props)
        self.assertTrue(out["drug_like"])
        self.assertEqual(out["violation_count"], 0)

    def test_violations_counted(self):
        props = {"molecular_weight": 800, "xlogp": 7,
                 "h_bond_donors": 8, "h_bond_acceptors": 15}
        out = pubchem._lipinski(props)
        self.assertEqual(out["violation_count"], 4)
        self.assertFalse(out["drug_like"])


class TestExtendedRules(unittest.TestCase):
    def test_none_when_descriptors_missing(self):
        out = pubchem._extended_rules({})
        self.assertIsNone(out["veber"]["pass"])
        self.assertIsNone(out["egan"]["pass"])

    def test_passing_small_molecule(self):
        props = {"molecular_weight": 200, "xlogp": 1.0,
                 "tpsa": 40, "rotatable_bonds": 3}
        out = pubchem._extended_rules(props)
        self.assertTrue(out["veber"]["pass"])
        self.assertTrue(out["lead_like"]["pass"])
        self.assertEqual(out["absorption"]["gi_absorption"], "High")


if __name__ == "__main__":
    unittest.main()
