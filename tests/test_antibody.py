"""Offline tests for the antibody layer (typing, CDRs, liabilities, interface)."""

import unittest

from snaclex import antibody, interactions
from snaclex.pdbparse import Atom, Structure
from tests.fixtures import structure as make_structure

# Real variable-domain sequences (trastuzumab) + a non-antibody control.
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
TRASTUZUMAB_VL = (
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFTGSR"
    "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
# Hen egg-white lysozyme — Ig-unrelated fold, has cysteines (false-positive guard).
LYSOZYME = (
    "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCND"
    "GRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"
)

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}


def _chain_atoms(seq, chain, start=1, x0=0.0):
    atoms = []
    for i, aa in enumerate(seq):
        atoms.append(
            Atom(
                serial=i + 1, name="CA", res_name=_ONE_TO_THREE[aa], chain=chain,
                res_seq=start + i, icode="", x=x0 + i * 3.8, y=0.0, z=0.0,
                element="C", is_hetero=False,
            )
        )
    return atoms


def _seq_structure(*chain_seqs):
    protein = []
    for chain, seq in chain_seqs:
        protein.extend(_chain_atoms(seq, chain))
    return make_structure(protein=protein)


class TestAntibodyTyping(unittest.TestCase):
    def test_detects_heavy_and_light_variable_chains(self):
        struct = _seq_structure(("H", TRASTUZUMAB_VH), ("L", TRASTUZUMAB_VL))
        ab = antibody.analyze(struct)
        self.assertTrue(ab["is_antibody"])
        types = {c["chain"]: c["type"] for c in ab["chains"]}
        self.assertEqual(types["H"], "VH")
        self.assertEqual(types["L"], "VL")

    def test_cdrs_delimited_for_variable_chains(self):
        struct = _seq_structure(("H", TRASTUZUMAB_VH))
        ab = antibody.analyze(struct)
        heavy = next(c for c in ab["chains"] if c["chain"] == "H")
        self.assertEqual(set(heavy["cdrs"]), {"CDR1", "CDR2", "CDR3"})
        cdr3 = heavy["cdrs"]["CDR3"]
        self.assertLess(cdr3["range"][0], cdr3["range"][1])
        self.assertTrue(ab["cdr_residues"])

    def test_non_antibody_not_flagged(self):
        struct = _seq_structure(("A", LYSOZYME))
        ab = antibody.analyze(struct)
        self.assertFalse(ab["is_antibody"])
        for c in ab["chains"]:
            self.assertNotIn(c["type"], ("VH", "VL"))


class TestAntibodyLiabilities(unittest.TestCase):
    def test_liability_scan_shape(self):
        struct = _seq_structure(("H", TRASTUZUMAB_VH))
        ab = antibody.analyze(struct)
        for l in ab["liabilities"]:
            self.assertIsInstance(l["in_cdr"], bool)
            self.assertIn("chain", l)
            self.assertIn("res_seq", l)
        self.assertEqual(ab["liability_count"], len(ab["liabilities"]))

    def test_deamidation_motif_detected(self):
        seq = "NG" + TRASTUZUMAB_VH  # explicit N-G hotspot at the start
        struct = _seq_structure(("H", seq))
        ab = antibody.analyze(struct)
        self.assertTrue(
            any(l["type"] == "Deamidation" and l["motif"] == "NG"
                for l in ab["liabilities"])
        )


class TestAntibodyInterface(unittest.TestCase):
    def test_interface_available_flag_with_antigen(self):
        struct = _seq_structure(
            ("H", TRASTUZUMAB_VH), ("L", TRASTUZUMAB_VL), ("A", LYSOZYME)
        )
        ab = antibody.analyze(struct)
        self.assertTrue(ab["interface_available"])

    def test_profile_interface_paratope_epitope_and_bsa(self):
        ab_atoms, ag_atoms = [], []
        for i in range(6):
            ab_atoms.append(Atom(i + 1, "CA", "ALA", "H", i + 1, "", i * 3.8, 0.0, 0.0, "C", False))
            ag_atoms.append(Atom(i + 1, "CA", "GLY", "A", i + 1, "", i * 3.8, 3.5, 0.0, "C", False))
        struct = Structure(
            atoms=ab_atoms + ag_atoms, protein_atoms=ab_atoms + ag_atoms,
            components=[], chains=["A", "H"],
        )
        result = interactions.profile_interface(struct, ["H"], ["A"])
        self.assertTrue(result["available"])
        self.assertGreater(result["paratope_residue_count"], 0)
        self.assertGreater(result["epitope_residue_count"], 0)
        self.assertIsNotNone(result["buried_surface_area_A2"])
        self.assertGreaterEqual(result["buried_surface_area_A2"], 0)


if __name__ == "__main__":
    unittest.main()
