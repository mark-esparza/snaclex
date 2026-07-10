"""Tests for the evidence model + A–F level invariants. Offline."""

import unittest

from snaclex import evidence
from snaclex.evidence import EvidenceInvariantError


class TestInvariants(unittest.TestCase):
    def test_docking_units_forced_to_score(self):
        ev = evidence.docking_evidence(
            protein_ref="P", chemical_ref="C", score=-7.2, software="SnaCleX",
            params={"seed": 42}, seed=42, structure_ref="1IEP")
        self.assertEqual(ev["level"], "E")
        self.assertEqual(ev["units"], "score")

    def test_docking_cannot_be_affinity(self):
        with self.assertRaises(EvidenceInvariantError):
            evidence.make_evidence(
                subject_ref="P", predicate="binds", object_ref="C", level="E",
                source="x", assay_type="Ki", computation={"software": "x"})

    def test_docking_requires_computation(self):
        with self.assertRaises(EvidenceInvariantError):
            evidence.make_evidence(
                subject_ref="P", predicate="binds", object_ref="C", level="E",
                source="x", value=-5, units="score")

    def test_homology_requires_transfer_basis(self):
        with self.assertRaises(EvidenceInvariantError):
            evidence.make_evidence(
                subject_ref="P", predicate="binds", object_ref="C", level="D",
                source="homology")

    def test_ml_requires_basis(self):
        with self.assertRaises(EvidenceInvariantError):
            evidence.make_evidence(
                subject_ref="P", predicate="binds", object_ref="C", level="F",
                source="model")

    def test_unknown_level_rejected(self):
        with self.assertRaises(EvidenceInvariantError):
            evidence.make_evidence(subject_ref="P", predicate="binds",
                                   object_ref="C", level="Z", source="x")


class TestBuilders(unittest.TestCase):
    def test_level_a_from_pdb(self):
        ev = evidence.level_a_from_pdb_ligand(
            protein_ref="P", chemical_ref="C", pdb_id="1IEP", chain="A",
            ligand_id="STI", method="X-RAY")
        self.assertEqual(ev["level"], "A")
        self.assertTrue(ev["is_direct"])
        self.assertIn("1IEP", ev["source_record"])

    def test_level_b_keeps_inactive(self):
        row = {"aid": 42, "outcome": "Inactive", "activity_name": "IC50",
               "activity_value_uM": None, "assay_type": "confirmatory", "pmid": None}
        ev = evidence.level_b_from_bioassay(
            protein_ref="P", chemical_ref="C", assay_row=row)
        self.assertEqual(ev["level"], "B")
        self.assertEqual(ev["outcome"], "Inactive")

    def test_level_d_carries_identity(self):
        ev = evidence.level_d_transfer(
            protein_ref="P", chemical_ref="C", from_protein="ABL1",
            sequence_identity=0.62, aligned_site_identity=0.9)
        self.assertEqual(ev["level"], "D")
        self.assertEqual(ev["transfer"]["sequence_identity"], 0.62)


class TestChemicalMatch(unittest.TestCase):
    def test_exact_cid(self):
        self.assertEqual(
            evidence.classify_chemical_match({"cid": 5291}, {"cid": 5291}), "exact")

    def test_stereo_same_skeleton(self):
        q = {"inchikey": "AAAAAAAAAAAAAA-BBBBBBBBBB-N"}
        c = {"inchikey": "AAAAAAAAAAAAAA-CCCCCCCCCC-N"}
        self.assertEqual(evidence.classify_chemical_match(q, c), "stereo")

    def test_parent(self):
        q = {"cid": 10, "parent_cid": 5}
        c = {"cid": 5}
        self.assertEqual(evidence.classify_chemical_match(q, c), "parent")

    def test_unresolved(self):
        self.assertEqual(
            evidence.classify_chemical_match({"cid": 1}, {"cid": 2}), "unresolved")


class TestSummaryNeverBoolean(unittest.TestCase):
    def test_counts_and_split(self):
        items = [
            evidence.level_a_from_pdb_ligand(
                protein_ref="P", chemical_ref="C", pdb_id="1IEP", chain="A",
                ligand_id="STI"),
            evidence.docking_evidence(
                protein_ref="P", chemical_ref="C", score=-5, software="SnaCleX",
                params={}, seed=1),
        ]
        s = evidence.summarize(items)
        self.assertEqual(s["by_level"]["A"], 1)
        self.assertEqual(s["by_level"]["E"], 1)
        self.assertEqual(s["direct_experimental"], 1)
        self.assertEqual(s["predicted_or_inferred"], 1)
        self.assertNotIn("interacts", s)  # never a boolean claim

    def test_rank_orders_by_directness(self):
        e_dock = evidence.docking_evidence(
            protein_ref="P", chemical_ref="C", score=-5, software="x", params={})
        e_pdb = evidence.level_a_from_pdb_ligand(
            protein_ref="P", chemical_ref="C", pdb_id="1IEP", chain="A", ligand_id="STI")
        ranked = evidence.rank([e_dock, e_pdb])
        self.assertEqual(ranked[0]["level"], "A")
        self.assertEqual(ranked[1]["level"], "E")


if __name__ == "__main__":
    unittest.main()
