"""Tests for protein-variant parsing and structure mapping (offline)."""

import unittest

from snaclex import variants
from tests.fixtures import atom, structure


def _chain(seq, *, chain="A", start=1):
    """Build a single-residue-per-position protein chain from a 1-letter seq."""
    one_to_three = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
        "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
        "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
    }
    atoms = []
    for i, aa in enumerate(seq):
        atoms.append(
            atom("C", float(i), 0.0, 0.0, name="CA",
                 res_name=one_to_three[aa], chain=chain, res_seq=start + i)
        )
    return atoms


class TestParseVariant(unittest.TestCase):
    def test_one_letter(self):
        v = variants.parse_variant("R273H")
        self.assertEqual((v["wt"], v["position"], v["mut"]), ("R", 273, "H"))

    def test_three_letter_with_prefix(self):
        v = variants.parse_variant("p.Arg273His")
        self.assertEqual((v["wt"], v["position"], v["mut"]), ("R", 273, "H"))

    def test_bare_position(self):
        v = variants.parse_variant("273")
        self.assertEqual(v["position"], 273)
        self.assertIsNone(v["wt"])

    def test_nonsense(self):
        v = variants.parse_variant("R273*")
        self.assertEqual(v["mut"], "*")

    def test_unparseable(self):
        self.assertIsNone(variants.parse_variant("not a variant"))
        self.assertIsNone(variants.parse_variant(""))


class TestMapping(unittest.TestCase):
    def setUp(self):
        # Structure sequence = "ACDEFGHIKL" at residues 1..10.
        self.seq = "ACDEFGHIKL"
        self.struct = structure(_chain(self.seq))

    def test_identity_mapping(self):
        # UniProt seq identical to structure -> position p maps to residue p.
        out = variants.annotate(self.struct, ["C2D", "K9A"], self.seq)
        self.assertEqual(out["mapped_count"], 2)
        by_pos = {r["position"]: r for r in out["variants"]}
        self.assertEqual(by_pos[2]["res_id"], "A/CYS2")
        self.assertTrue(by_pos[2]["wt_matches_structure"])
        self.assertEqual(by_pos[9]["res_id"], "A/LYS9")

    def test_offset_uniprot(self):
        # UniProt has a 5-residue N-terminal extension the structure lacks; the
        # alignment should still place a UniProt position on the right residue.
        uniprot = "MMMMM" + self.seq  # positions 1..5 absent from structure
        out = variants.annotate(self.struct, ["D8E"], uniprot)
        r = out["variants"][0]
        # UniProt position 8 == structure 'D' (3rd structure residue).
        self.assertTrue(r["mapped"])
        self.assertEqual(r["res_name"], "ASP")
        self.assertTrue(r["wt_matches_structure"])

    def test_wt_mismatch_flagged(self):
        out = variants.annotate(self.struct, ["A2V"], self.seq)
        r = out["variants"][0]
        self.assertTrue(r["mapped"])
        # Structure residue 2 is CYS, not ALA -> WT does not match.
        self.assertFalse(r["wt_matches_structure"])

    def test_unmapped_position(self):
        out = variants.annotate(self.struct, ["R999H"], self.seq)
        r = out["variants"][0]
        self.assertFalse(r["mapped"])
        self.assertIn("not covered", r["reason"])

    def test_fasta_parse(self):
        text = ">sp|P04637|P53_HUMAN\nACDEF\nGHIKL\n"
        self.assertEqual(variants._parse_fasta(text), "ACDEFGHIKL")


if __name__ == "__main__":
    unittest.main()
