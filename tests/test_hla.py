"""Offline tests for the HLA / MHC module."""

import unittest

from snaclex import hla
from snaclex.pdbparse import Atom
from tests.fixtures import structure as make_structure

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}

LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCND"
    "GRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
)


def _chain_atoms(seq, chain, start=1):
    return [
        Atom(i + 1, "CA", _ONE_TO_THREE.get(a, "ALA"), chain, start + i, "",
             i * 3.8, 0, 0, "C", False)
        for i, a in enumerate(seq)
    ]


def _structure(*chain_seqs):
    protein = []
    for chain, seq in chain_seqs:
        protein.extend(_chain_atoms(seq, chain))
    return make_structure(protein=protein)


class TestAlleleNomenclature(unittest.TestCase):
    def test_normalize(self):
        for raw, want in [
            ("HLA-B*57:01", "B*57:01"),
            ("B*57:01", "B*57:01"),
            ("HLA-B5701", "B*57:01"),
            ("A*02:01", "A*02:01"),
            ("DRB1*15:01", "DRB1*15:01"),
        ]:
            self.assertEqual(hla.normalize_allele(raw), want)
        self.assertIsNone(hla.normalize_allele("not-an-allele"))

    def test_class(self):
        self.assertEqual(hla.allele_class("B*57:01"), "I")
        self.assertEqual(hla.allele_class("DRB1*15:01"), "II")

    def test_drug_hits(self):
        hits = hla.drug_hits("B*57:01")
        self.assertTrue(any(h["drug"] == "abacavir" for h in hits))
        self.assertEqual(hla.drug_hits("A*02:01"), [])
        # Every seeded association carries a citation.
        for allele, rows in hla.DRUG_ASSOCIATIONS.items():
            for r in rows:
                self.assertTrue(r.get("citation"))


class TestDetection(unittest.TestCase):
    def test_detect_class_i(self):
        struct = _structure(("A", hla.CLASS_I_REF), ("B", hla.B2M_REF), ("C", "SIINFEKL"))
        det = hla.detect(struct)
        self.assertTrue(det["is_mhc"])
        self.assertEqual(det["mhc_class"], "I")
        self.assertEqual(det["heavy_chain"], "A")
        self.assertEqual(det["b2m_chain"], "B")

    def test_non_mhc_not_detected(self):
        det = hla.detect(_structure(("A", LYSOZYME)))
        self.assertFalse(det["is_mhc"])


class TestGrooveAnalysis(unittest.TestCase):
    def test_groove_reference_has_no_differences(self):
        # Structure IS the reference allele -> zero groove differences.
        struct = _structure(("A", hla.CLASS_I_REF), ("B", hla.B2M_REF))
        result = hla.analyze(struct, "A*02:01")
        self.assertTrue(result["is_mhc"])
        self.assertIn("groove", result)
        self.assertEqual(result["groove"]["differences"], [])
        self.assertTrue(result["groove"]["groove_residues"])

    def test_groove_polymorphism_detected(self):
        # Mutate groove position 9 (pocket B/C) from F to something else.
        seq = list(hla.CLASS_I_REF)
        self.assertEqual(seq[8], "F")  # position 9 is Phe in A*02:01
        seq[8] = "D"  # F9D — a charge-introducing change
        struct = _structure(("A", "".join(seq)), ("B", hla.B2M_REF))
        result = hla.analyze(struct, "A*02:01")
        diffs = result["groove"]["differences"]
        self.assertTrue(any(d["res_seq"] == 9 and d["allele"] == "ASP" for d in diffs))

    def test_analyze_reports_drug_hits(self):
        struct = _structure(("A", hla.CLASS_I_REF), ("B", hla.B2M_REF))
        result = hla.analyze(struct, "HLA-B*57:01")
        self.assertEqual(result["allele"], "B*57:01")
        self.assertTrue(any(h["drug"] == "abacavir" for h in result["drug_associations"]))

    def test_peptide_docking_deferred(self):
        struct = _structure(("A", hla.CLASS_I_REF), ("B", hla.B2M_REF))
        result = hla.analyze(struct)
        self.assertFalse(result["peptide_docking_available"])


if __name__ == "__main__":
    unittest.main()
