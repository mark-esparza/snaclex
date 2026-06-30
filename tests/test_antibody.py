from __future__ import annotations

import unittest

from snaclex.antibody import detect_cdrs
from snaclex.pdbparse import Structure

from tests.fixtures import atom

AA1TO3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}

# Trastuzumab (Herceptin) VH domain. Kabat-annotated CDRs (for reference):
#   CDR-H1 ~ GFNIKDTYIH, CDR-H2 ~ RIYPTNGYTRYADSVKG, CDR-H3 ~ WGGDGFYAMDY
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRY"
    "ADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)


def _chain_structure(seq: str, chain: str = "H") -> Structure:
    atoms = [
        atom("C", float(i), 0.0, 0.0, name="CA", res_name=AA1TO3[aa],
             chain=chain, res_seq=i + 1, serial=i + 1)
        for i, aa in enumerate(seq)
    ]
    protein_atoms = list(atoms)
    return Structure(
        atoms=atoms, protein_atoms=protein_atoms, components=[], chains=[chain]
    )


class TestDetectCdrs(unittest.TestCase):
    def test_detects_three_cdrs_in_antibody_vh_domain(self):
        structure = _chain_structure(TRASTUZUMAB_VH)
        loops = detect_cdrs(structure)
        names = [loop.name for loop in loops]
        self.assertEqual(names, ["CDR1", "CDR2", "CDR3"])
        self.assertTrue(all(loop.chain == "H" for loop in loops))

    def test_cdr3_is_high_confidence_and_matches_known_loop(self):
        structure = _chain_structure(TRASTUZUMAB_VH)
        loops = {loop.name: loop for loop in detect_cdrs(structure)}
        cdr3 = loops["CDR3"]
        self.assertEqual(cdr3.confidence, "high")
        # The known CDR-H3 loop should fall inside the detected window.
        window = TRASTUZUMAB_VH[cdr3.start_res_seq - 1 : cdr3.end_res_seq]
        self.assertIn("WGGDGFYAMDY", window)

    def test_loops_are_ordered_and_non_overlapping(self):
        structure = _chain_structure(TRASTUZUMAB_VH)
        loops = detect_cdrs(structure)
        for a, b in zip(loops, loops[1:]):
            self.assertLess(a.end_res_seq, b.start_res_seq)

    def test_non_antibody_chain_yields_no_cdrs(self):
        # A short, non-Ig-sized chain shouldn't trigger detection.
        structure = _chain_structure("AGVL" * 5)
        self.assertEqual(detect_cdrs(structure), [])

    def test_empty_structure_yields_no_cdrs(self):
        structure = Structure(atoms=[], protein_atoms=[], components=[], chains=[])
        self.assertEqual(detect_cdrs(structure), [])


if __name__ == "__main__":
    unittest.main()
