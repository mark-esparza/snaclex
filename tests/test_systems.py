"""Tests for curated system definitions and their verification state."""

import json
import os
import tempfile
import unittest
from unittest import mock

from snaclex import systems


MINIMAL = {
    "id": "demo",
    "name": "Demo system",
    "targets": [
        {"id": "A", "pdb": "1ABC", "uniprot": "P00001", "site": {"ligand": "LIG"}},
        {"id": "B", "pdb": "2DEF", "site": {"pocket": 0}},
    ],
    "ligands": [
        {"id": "x", "query": "dopamine"},
        {"id": "y", "query": "epinephrine"},
    ],
}


def _doc(**overrides):
    doc = json.loads(json.dumps(MINIMAL))
    doc.update(overrides)
    return doc


class TestValidation(unittest.TestCase):
    def test_minimal_document_validates(self):
        doc = _doc()
        self.assertIs(systems.validate_system(doc), doc)

    def test_missing_targets_rejected(self):
        with self.assertRaises(systems.SystemError_):
            systems.validate_system(_doc(targets=[]))

    def test_missing_ligands_rejected(self):
        with self.assertRaises(systems.SystemError_):
            systems.validate_system(_doc(ligands=[]))

    def test_target_without_site_rejected(self):
        bad = _doc()
        bad["targets"][0].pop("site")
        with self.assertRaises(systems.SystemError_) as ctx:
            systems.validate_system(bad)
        self.assertIn("site", str(ctx.exception))

    def test_target_with_unknown_site_key_rejected(self):
        bad = _doc()
        bad["targets"][0]["site"] = {"magic": 1}
        with self.assertRaises(systems.SystemError_):
            systems.validate_system(bad)

    def test_duplicate_target_ids_rejected(self):
        bad = _doc()
        bad["targets"][1]["id"] = "A"
        with self.assertRaises(systems.SystemError_) as ctx:
            systems.validate_system(bad)
        self.assertIn("duplicate", str(ctx.exception))

    def test_duplicate_ligand_ids_rejected(self):
        bad = _doc()
        bad["ligands"][1]["id"] = "x"
        with self.assertRaises(systems.SystemError_):
            systems.validate_system(bad)

    def test_target_missing_pdb_rejected(self):
        bad = _doc()
        bad["targets"][0].pop("pdb")
        with self.assertRaises(systems.SystemError_):
            systems.validate_system(bad)


class TestVerificationSummary(unittest.TestCase):
    def test_unverified_system_is_labelled_a_draft(self):
        summary = systems.verification_summary(_doc())
        self.assertFalse(summary["fully_verified"])
        self.assertEqual(summary["n_verified"], 0)
        self.assertEqual(sorted(summary["unverified_target_ids"]), ["A", "B"])
        self.assertIn("DRAFT", summary["note"])

    def test_partially_verified_still_a_draft(self):
        doc = _doc()
        doc["targets"][0]["verified"] = True
        summary = systems.verification_summary(doc)
        self.assertFalse(summary["fully_verified"])
        self.assertEqual(summary["unverified_target_ids"], ["B"])

    def test_fully_verified_system_says_so(self):
        doc = _doc()
        for target in doc["targets"]:
            target["verified"] = True
        summary = systems.verification_summary(doc)
        self.assertTrue(summary["fully_verified"])
        self.assertNotIn("DRAFT", summary["note"])


class TestLoading(unittest.TestCase):
    def test_unknown_system_raises(self):
        with self.assertRaises(systems.SystemError_):
            systems.load_system("no-such-system")

    def test_path_traversal_rejected(self):
        for bad in ("../secrets", "a/b", "..", "x.json", ""):
            with self.assertRaises(systems.SystemError_):
                systems.load_system(bad)

    def test_shipped_catecholamine_system_loads_and_is_well_formed(self):
        doc = systems.load_system("catecholamine")
        self.assertEqual(doc["id"], "catecholamine")
        self.assertTrue(doc["targets"])
        self.assertTrue(doc["ligands"])
        # Within the engine's panel limits, so it can actually be run.
        from snaclex import panel
        self.assertLessEqual(len(doc["targets"]), panel.MAX_TARGETS)
        self.assertLessEqual(len(doc["ligands"]), panel.MAX_LIGANDS)
        self.assertLessEqual(
            len(doc["targets"]) * len(doc["ligands"]), panel.MAX_CELLS
        )

    def test_shipped_system_is_honest_about_being_unverified(self):
        """It ships as a draft; nothing should claim otherwise."""
        doc = systems.load_system("catecholamine")
        self.assertFalse(doc["verification"]["fully_verified"])
        self.assertIn("DRAFT", doc["verification"]["note"])
        for target in doc["targets"]:
            self.assertIsNot(target.get("verified"), True)

    def test_catalog_lists_the_shipped_system(self):
        ids = {entry["id"] for entry in systems.list_systems()}
        self.assertIn("catecholamine", ids)

    def test_catalog_entries_carry_verification_state(self):
        for entry in systems.list_systems():
            self.assertIn("verification", entry)
            self.assertIn("fully_verified", entry["verification"])


class TestPanelSpecs(unittest.TestCase):
    def test_unfiltered_returns_everything(self):
        targets, ligands = systems.panel_specs(_doc())
        self.assertEqual(len(targets), 2)
        self.assertEqual(len(ligands), 2)

    def test_filters_select_a_slice(self):
        targets, ligands = systems.panel_specs(
            _doc(), target_ids=["B"], ligand_ids=["y"]
        )
        self.assertEqual([t["id"] for t in targets], ["B"])
        self.assertEqual([lig["id"] for lig in ligands], ["y"])

    def test_filter_matching_nothing_raises(self):
        with self.assertRaises(systems.SystemError_):
            systems.panel_specs(_doc(), target_ids=["nope"])
        with self.assertRaises(systems.SystemError_):
            systems.panel_specs(_doc(), ligand_ids=["nope"])


class TestVerify(unittest.TestCase):
    """verify_system with injected fakes — no network."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = os.path.join(self.tmp.name, "demo.json")
        with open(self.path, "w", encoding="utf-8") as fh:
            json.dump(_doc(), fh)
        patcher = mock.patch.object(systems, "_DATA_DIR", self.tmp.name)
        patcher.start()
        self.addCleanup(patcher.stop)

    def _run(self, *, ligands_present=("LIG",), accessions=("P00001",), **kw):
        from snaclex.pdbparse import Component

        def fake_meta(pdb):
            return {"title": f"{pdb} title", "resolution_A": 1.9,
                    "experimental_method": "X-RAY DIFFRACTION"}

        def fake_struct(pdb):
            return "TEXT"

        def fake_parse(text):
            comps = [Component(name, "A", 1, "", [])
                     for name in ligands_present]
            extra = [Component(name, "A", 2, "", [])
                     for name in kw.pop("cofactors_present", ())]
            return mock.Mock(ligand_components=comps, components=comps + extra)

        return systems.verify_system(
            "demo",
            fetch_metadata=fake_meta,
            fetch_structure=fake_struct,
            parse_structure=fake_parse,
            fetch_uniprots=lambda pdb: list(accessions),
            **kw,
        )

    def test_clean_verification_marks_targets_verified(self):
        result = self._run()
        self.assertEqual(result["n_problem"], 0)
        self.assertEqual(result["n_ok"], 2)

    def test_missing_site_ligand_is_reported_as_a_problem(self):
        result = self._run(ligands_present=("XXX",))
        problems = [f for f in result["findings"] if f["problems"]]
        self.assertEqual(len(problems), 1)
        self.assertIn("not found", problems[0]["problems"][0])

    def test_declared_cofactor_present_is_recorded(self):
        doc = _doc()
        doc["targets"][0]["cofactors"] = ["FAD"]
        with open(self.path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh)
        result = self._run(ligands_present=("LIG", "FAD"))
        self.assertEqual(result["findings"][0]["cofactors_found"], ["FAD"])
        self.assertEqual(result["n_problem"], 0)

    def test_missing_cofactor_is_reported_as_a_problem(self):
        doc = _doc()
        doc["targets"][0]["cofactors"] = ["FAD"]
        with open(self.path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh)
        result = self._run(ligands_present=("LIG",))
        problems = result["findings"][0]["problems"]
        self.assertTrue(problems)
        self.assertIn("empty cavity", problems[0])

    def test_any_one_declared_cofactor_alternative_satisfies(self):
        """SAM or SAH both count: entries differ in which is deposited."""
        doc = _doc()
        doc["targets"][0]["cofactors"] = ["SAM", "SAH"]
        with open(self.path, "w", encoding="utf-8") as fh:
            json.dump(doc, fh)
        result = self._run(ligands_present=("LIG", "SAH"))
        self.assertEqual(result["findings"][0]["cofactors_found"], ["SAH"])
        self.assertEqual(result["n_problem"], 0)

    def test_wrong_uniprot_is_reported_as_a_problem(self):
        result = self._run(accessions=("Q99999",))
        problems = [f for f in result["findings"] if f["problems"]]
        self.assertTrue(problems)
        self.assertIn("UniProt", problems[0]["problems"][0])

    def test_fetch_failure_is_reported_not_raised(self):
        def boom(pdb):
            raise RuntimeError("offline")

        result = systems.verify_system(
            "demo", fetch_metadata=boom,
            fetch_structure=lambda p: "", parse_structure=lambda t: None,
            fetch_uniprots=lambda p: [],
        )
        self.assertEqual(result["n_problem"], 2)
        self.assertIn("could not fetch entry", result["findings"][0]["problems"][0])

    def _read(self):
        with open(self.path, encoding="utf-8") as fh:
            return fh.read()

    def test_dry_run_does_not_touch_the_file(self):
        before = self._read()
        self._run()
        self.assertEqual(self._read(), before)

    def test_write_persists_verified_flags(self):
        self._run(write=True)
        saved = json.loads(self._read())
        self.assertTrue(all(t["verified"] for t in saved["targets"]))
        self.assertTrue(saved["last_verified_utc"])
        # The derived block is never written back into the source of truth.
        self.assertNotIn("verification", saved)

    def test_written_file_still_loads_and_reports_verified(self):
        self._run(write=True)
        doc = systems.load_system("demo")
        self.assertTrue(doc["verification"]["fully_verified"])


if __name__ == "__main__":
    unittest.main()
