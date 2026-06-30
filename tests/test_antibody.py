"""Tests for snaclex.antibody — CDR loop annotation (offline)."""

import unittest
from types import SimpleNamespace

from snaclex import antibody


def _make_residue(chain, res_seq, total=3, types=None, min_distance=3.2):
    return {
        "res_id": f"{chain}/{res_seq}",
        "chain": chain,
        "res_seq": res_seq,
        "res_name": "ALA",
        "types": types or ["hydrogen_bond"],
        "total": total,
        "min_distance": min_distance,
    }


def _make_profile(res_a, res_b=None):
    return {
        "mode": "protein-protein",
        "chains_a": ["H", "L"],
        "chains_b": ["A"],
        "interactions": [],
        "counts": {"hydrogen_bond": 4, "salt_bridge": 1, "hydrophobic": 2, "aromatic": 0},
        "interaction_total": 7,
        "interface_residues_a": res_a,
        "interface_residues_b": res_b or [],
        "interface_residue_count": len(res_a) + len(res_b or []),
    }


class TestCdrBoundaries(unittest.TestCase):
    def test_kabat_heavy_h1(self):
        res = [_make_residue("H", 32)]  # inside H1 31-35
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "H1")

    def test_kabat_heavy_h2(self):
        res = [_make_residue("H", 55)]  # inside H2 50-65
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "H2")

    def test_kabat_heavy_h3(self):
        res = [_make_residue("H", 99)]  # inside H3 95-102
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "H3")

    def test_kabat_heavy_framework(self):
        res = [_make_residue("H", 40)]  # between H1 and H2 — framework
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertIsNone(p["interface_residues_a"][0]["cdr"])

    def test_kabat_light_l1(self):
        res = [_make_residue("L", 28)]  # inside L1 24-34
        p = antibody.annotate_paratope(_make_profile(res), heavy=None, light="L")
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "L1")

    def test_kabat_light_l2(self):
        res = [_make_residue("L", 52)]  # inside L2 50-56
        p = antibody.annotate_paratope(_make_profile(res), heavy=None, light="L")
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "L2")

    def test_kabat_light_l3(self):
        res = [_make_residue("L", 93)]  # inside L3 89-97
        p = antibody.annotate_paratope(_make_profile(res), heavy=None, light="L")
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "L3")

    def test_imgt_scheme(self):
        res = [_make_residue("H", 30)]  # inside IMGT H1 27-38
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None, scheme="imgt")
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "H1")
        self.assertEqual(p["antibody_chains"]["scheme"], "imgt")

    def test_chothia_h1(self):
        res = [_make_residue("H", 29)]  # inside Chothia H1 26-32
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None, scheme="chothia")
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "H1")

    def test_unknown_scheme_falls_back_to_kabat(self):
        res = [_make_residue("H", 32)]
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None, scheme="bogus")
        self.assertEqual(p["antibody_chains"]["scheme"], "kabat")
        self.assertEqual(p["interface_residues_a"][0]["cdr"], "H1")

    def test_chain_not_heavy_or_light_gives_none(self):
        res = [_make_residue("X", 32)]
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light="L")
        self.assertIsNone(p["interface_residues_a"][0]["cdr"])


class TestCdrSummary(unittest.TestCase):
    def test_summary_counts(self):
        res = [
            _make_residue("H", 32, total=2),  # H1
            _make_residue("H", 33, total=3),  # H1
            _make_residue("H", 99, total=1),  # H3
        ]
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertEqual(p["cdr_summary"]["H1"], 5)
        self.assertEqual(p["cdr_summary"]["H3"], 1)
        self.assertNotIn("H2", p["cdr_summary"])

    def test_has_cdr_annotation_flag(self):
        p = antibody.annotate_paratope(_make_profile([]), heavy="H", light=None)
        self.assertTrue(p["has_cdr_annotation"])

    def test_antibody_chains_recorded(self):
        p = antibody.annotate_paratope(_make_profile([]), heavy="H", light="L", scheme="imgt")
        self.assertEqual(p["antibody_chains"], {"heavy": "H", "light": "L", "scheme": "imgt"})

    def test_original_profile_not_mutated(self):
        res = [_make_residue("H", 32)]
        original = _make_profile(res)
        original_a = list(original["interface_residues_a"])
        antibody.annotate_paratope(original, heavy="H", light=None)
        self.assertNotIn("cdr", original["interface_residues_a"][0])
        self.assertEqual(original["interface_residues_a"], original_a)


class TestNumberingNote(unittest.TestCase):
    def test_warns_on_high_res_seq(self):
        res = [_make_residue("H", 250)]  # far beyond Kabat range
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertIsNotNone(p["numbering_note"])
        self.assertIn("250", p["numbering_note"])

    def test_no_warning_within_kabat_range(self):
        res = [_make_residue("H", 99)]
        p = antibody.annotate_paratope(_make_profile(res), heavy="H", light=None)
        self.assertIsNone(p["numbering_note"])

    def test_no_warning_when_no_chains_given(self):
        res = [_make_residue("H", 99)]
        p = antibody.annotate_paratope(_make_profile(res), heavy=None, light=None)
        self.assertIsNone(p["numbering_note"])


class TestDetectVhvl(unittest.TestCase):
    def _make_structure(self, chain_res_seqs: dict):
        """Build a minimal structure mock with protein_atoms."""
        atoms = []
        for chain, seqs in chain_res_seqs.items():
            for seq in seqs:
                atoms.append(SimpleNamespace(chain=chain, res_seq=seq))
        return SimpleNamespace(protein_atoms=atoms)

    def test_fab_length_detection(self):
        # VH+CH1 ~220 residues, VL+CL ~210 residues (Kabat-numbered)
        s = self._make_structure({
            "H": list(range(1, 221)),   # 220 residues
            "L": list(range(1, 211)),   # 210 residues
            "A": list(range(1, 400)),   # antigen
        })
        det = antibody.detect_vhvl_chains(s)
        self.assertEqual(det["confidence"], "high")
        self.assertEqual(det["heavy"], "H")
        self.assertEqual(det["light"], "L")

    def test_vd_only_fallback(self):
        # Two variable-domain-length chains (~120 residues each)
        s = self._make_structure({
            "H": list(range(1, 121)),
            "L": list(range(1, 121)),
        })
        det = antibody.detect_vhvl_chains(s)
        self.assertEqual(det["confidence"], "low")
        # Both detected but tentative
        self.assertIsNotNone(det["heavy"])
        self.assertIsNotNone(det["light"])

    def test_no_antibody_chains(self):
        # A random 300-residue chain
        s = self._make_structure({"A": list(range(1, 301))})
        det = antibody.detect_vhvl_chains(s)
        self.assertIsNone(det["heavy"])
        self.assertIsNone(det["light"])
        self.assertEqual(det["confidence"], "low")

    def test_picks_kabat_covered_chain(self):
        # Two Fab-length chains; one covers Kabat CDR windows, one uses sequential numbers
        # H: 1..220 (Kabat-like, CDR windows covered)
        # H2: 221..440 (sequential, CDR windows not covered)
        s = self._make_structure({
            "H": list(range(1, 221)),
            "H2": list(range(221, 441)),
            "L": list(range(1, 211)),
        })
        det = antibody.detect_vhvl_chains(s)
        self.assertEqual(det["heavy"], "H")


if __name__ == "__main__":
    unittest.main()
