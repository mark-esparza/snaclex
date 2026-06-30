"""Tests for HLA detection, groove derivation, and variant classification.

Builds a synthetic class-I-like complex offline: a long heavy chain (A), a
beta-2-microglobulin-sized chain (B), and a 9-residue peptide (C) positioned so
that specific heavy-chain residues contact peptide P2 and the C-terminus (the
anchor pockets) plus one mid-peptide position (groove-lining only).
"""

import unittest

from snaclex import hla
from tests.fixtures import atom, structure

# Heavy-chain residues we place against the peptide (the rest sit far away).
_GROOVE = {100, 120, 150}


def _hla_structure():
    protein = []
    # Heavy chain A: 270 residues; all far away except the groove residues.
    for i in range(1, 271):
        if i in _GROOVE:
            continue
        protein.append(atom("C", 500.0 + i, 0.0, 0.0, name="CA",
                             res_name="ALA", chain="A", res_seq=i))
    # Groove heavy residues near the peptide.
    protein.append(atom("N", 3.0, 0.0, 0.0, name="N", res_name="GLN",
                        chain="A", res_seq=100))   # contacts peptide P2
    protein.append(atom("O", 3.0, 10.0, 0.0, name="O", res_name="SER",
                        chain="A", res_seq=120))   # contacts peptide P5 (mid)
    protein.append(atom("O", 3.0, 20.0, 0.0, name="OD1", res_name="ASP",
                        chain="A", res_seq=150))   # contacts peptide P9 (C-term)
    # beta-2-microglobulin chain B: 99 residues, far away.
    for i in range(1, 100):
        protein.append(atom("C", -500.0 - i, 0.0, 0.0, name="CA",
                             res_name="ALA", chain="B", res_seq=i))
    # Peptide chain C: 9 residues 1..9; backbone O placed for 3 of them.
    pep_coords = {2: (0.0, 0.0, 0.0), 5: (0.0, 10.0, 0.0), 9: (0.0, 20.0, 0.0)}
    for i in range(1, 10):
        x, y, z = pep_coords.get(i, (100.0, 100.0 + i, 0.0))
        protein.append(atom("O", x, y, z, name="O", res_name="GLY",
                             chain="C", res_seq=i))
    return structure(protein)


class TestDetect(unittest.TestCase):
    def test_detects_class_i(self):
        det = hla.detect(_hla_structure())
        self.assertTrue(det["is_hla"])
        self.assertEqual(det["mhc_class"], "I")
        self.assertEqual(det["heavy_chain"], "A")
        self.assertEqual(det["b2m_chain"], "B")
        self.assertEqual(det["peptide_chain"], "C")

    def test_title_signal_raises_confidence(self):
        det = hla.detect(_hla_structure(), title="Crystal structure of HLA-A*02:01")
        self.assertEqual(det["confidence"], "high")
        self.assertTrue(any("title" in s for s in det["signals"]))

    def test_non_hla(self):
        prot = [atom("C", float(i), 0, 0, name="CA", res_name="ALA",
                     chain="A", res_seq=i) for i in range(1, 30)]
        det = hla.detect(structure(prot))
        self.assertFalse(det["is_hla"])


class TestGroove(unittest.TestCase):
    def setUp(self):
        self.struct = _hla_structure()
        self.det = hla.detect(self.struct)
        self.groove = hla.analyze_groove(self.struct, self.det)

    def test_groove_available(self):
        self.assertTrue(self.groove["available"])
        self.assertEqual(self.groove["peptide_length"], 9)

    def test_anchor_positions_are_p2_and_cterm(self):
        self.assertEqual(self.groove["anchor_peptide_positions"], [2, 9])

    def test_groove_and_anchor_residues(self):
        ids = set(self.groove["groove_residue_ids"])
        self.assertEqual(ids, {"A/GLN100", "A/SER120", "A/ASP150"})
        # Anchor-pocket residues are those contacting peptide P2 and P9.
        self.assertEqual(
            set(self.groove["anchor_pocket_residue_ids"]), {"A/GLN100", "A/ASP150"}
        )


class TestClassifyVariants(unittest.TestCase):
    def test_classification(self):
        struct = _hla_structure()
        det = hla.detect(struct)
        groove = hla.analyze_groove(struct, det)
        mapped = [
            {"input": "X100Y", "res_id": "A/GLN100", "mapped": True, "wt": "X", "mut": "Y"},
            {"input": "X120Y", "res_id": "A/SER120", "mapped": True, "wt": "X", "mut": "Y"},
            {"input": "X9Y", "res_id": "A/ALA9", "mapped": True, "wt": "X", "mut": "Y"},
            {"input": "bad", "mapped": False},
        ]
        out = hla.classify_variants(groove, mapped)
        by_res = {r["res_id"]: r["groove_class"] for r in out}
        self.assertEqual(by_res["A/GLN100"], "anchor-pocket")
        self.assertEqual(by_res["A/SER120"], "groove-lining")
        self.assertEqual(by_res["A/ALA9"], "peripheral")
        # Unmapped variants are dropped.
        self.assertEqual(len(out), 3)


if __name__ == "__main__":
    unittest.main()
