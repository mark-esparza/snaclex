"""Tests for the many-to-many ligand x target panel engine.

All I/O is injected, so these run fully offline: fake loaders hand the engine
small synthetic receptors and ligands, and the real docking/interaction code
does the rest.
"""

import unittest

from snaclex import docking, panel
from tests.fixtures import atom, structure


def _receptor(offset=0.0):
    """A small protein-atom cluster; `offset` shifts the whole site in x."""
    coords = [
        ("C", 3, 0, 0), ("O", -3, 0, 0), ("N", 0, 3, 0), ("C", 0, -3, 0),
        ("C", 2, 2, 1), ("O", -2, 2, -1), ("N", 2, -2, 1), ("C", -2, -2, -1),
        ("C", 4, 1, 0), ("C", -4, -1, 0), ("O", 1, 4, 0), ("N", -1, -4, 0),
    ]
    return structure(
        [atom(el, x + offset, y, z, name=el, res_name="LEU", res_seq=i)
         for i, (el, x, y, z) in enumerate(coords)]
    )


def _ligand(n_atoms=3, element="C"):
    return [
        {"element": element, "x": 1.4 * i, "y": 0.0, "z": 0.0}
        for i in range(n_atoms)
    ]


class PanelHarness:
    """Injectable fakes recording what the engine asked for."""

    def __init__(self, targets=2, fail_ligand=None, fail_target=None):
        self.structures = {
            f"T{i}": _receptor(offset=i * 0.5) for i in range(targets)
        }
        self.fail_ligand = fail_ligand
        self.fail_target = fail_target
        self.grid_builds = []
        self.ligand_loads = []

    def load_structure(self, pdb):
        if pdb == self.fail_target:
            raise RuntimeError("structure unavailable")
        return ("TEXT", self.structures[pdb], {"title": f"{pdb} title"})

    def load_ligand(self, query):
        self.ligand_loads.append(query)
        if query == self.fail_ligand:
            raise RuntimeError("no such compound")
        size = {"small": 2, "medium": 4, "large": 6}.get(query, 3)
        return {
            "atoms": _ligand(size),
            "source": "3d",
            "cid": 100 + size,
            "name": query,
            "formula": "C6H6",
        }

    def resolve_site(self, pdb, struct, target):
        return ((0.0, 0.0, 0.0), f"{pdb} test site")

    def build_grid(self, pdb, struct, center, extra_atoms=None):
        self.grid_builds.append((pdb, center, list(extra_atoms or [])))
        return docking.build_grid(struct, center, extra_atoms)

    def run(self, targets, ligands, **kw):
        return panel.run_panel(
            targets, ligands,
            load_structure=self.load_structure,
            load_ligand=self.load_ligand,
            resolve_site=self.resolve_site,
            build_grid=self.build_grid,
            seeds=kw.pop("seeds", 8),
            mc_steps=kw.pop("mc_steps", 5),
            **kw,
        )


def _targets(n):
    return [{"id": f"T{i}", "pdb": f"T{i}", "label": f"Target {i}"}
            for i in range(n)]


def _ligands(names):
    return [{"id": n, "query": n, "label": n} for n in names]


class TestPanelShape(unittest.TestCase):
    def test_every_pair_produces_a_cell(self):
        h = PanelHarness(targets=2)
        out = h.run(_targets(2), _ligands(["small", "medium", "large"]))
        self.assertEqual(out["n_cells"], 6)
        self.assertEqual(out["n_scored"], 6)
        self.assertEqual(out["n_failed"], 0)
        pairs = {(c["target_id"], c["ligand_id"]) for c in out["cells"]}
        self.assertEqual(len(pairs), 6)

    def test_grid_built_once_per_target_not_per_cell(self):
        """The whole reason M x N is affordable: N grids, not M*N."""
        h = PanelHarness(targets=2)
        h.run(_targets(2), _ligands(["small", "medium", "large"]))
        self.assertEqual(len(h.grid_builds), 2)

    def test_ligand_fetched_once_per_panel(self):
        h = PanelHarness(targets=3)
        h.run(_targets(3), _ligands(["small", "medium"]))
        self.assertEqual(sorted(h.ligand_loads), ["medium", "small"])

    def test_private_handles_not_leaked_into_result(self):
        """Result must be JSON-serializable; Grid/Structure objects are not."""
        import json

        h = PanelHarness(targets=2)
        out = h.run(_targets(2), _ligands(["small", "medium"]))
        for entry in out["targets"]:
            self.assertNotIn("_grid", entry)
            self.assertNotIn("_structure", entry)
        for entry in out["ligands"]:
            self.assertNotIn("_atoms", entry)
        json.dumps(out)  # raises if anything non-serializable survived


class TestSiteSubset(unittest.TestCase):
    """Trimming the receptor to the site shell must not change any result."""

    def test_distant_residues_are_dropped(self):
        near = [atom("C", 1, 0, 0, res_name="LEU", chain="A", res_seq=1)]
        far = [atom("C", 500, 0, 0, res_name="LEU", chain="A", res_seq=2)]
        s = structure(near + far)
        sub = panel.site_subset(s, (0, 0, 0))
        self.assertEqual(len(sub.protein_atoms), 1)
        self.assertEqual(sub.protein_atoms[0].res_seq, 1)

    def test_residues_are_kept_whole(self):
        """A residue straddling the radius is kept entire, so rings survive."""
        radius = panel.SITE_RADIUS
        res = [
            atom("C", radius - 1, 0, 0, res_name="PHE", chain="A", res_seq=7),
            atom("C", radius + 5, 0, 0, res_name="PHE", chain="A", res_seq=7),
        ]
        sub = panel.site_subset(structure(res), (0, 0, 0))
        self.assertEqual(len(sub.protein_atoms), 2)

    def test_empty_shell_falls_back_to_the_full_structure(self):
        s = structure([atom("C", 999, 0, 0, res_name="LEU", chain="A", res_seq=1)])
        self.assertIs(panel.site_subset(s, (0, 0, 0)), s)

    def test_radius_exceeds_box_reach_plus_widest_cutoff(self):
        from snaclex import interactions as inter
        widest = max(inter.HB_MAX, inter.SALT_MAX, inter.HYDRO_MAX,
                     inter.METAL_MAX, inter.ARO_CENTROID_MAX)
        self.assertGreater(panel.SITE_RADIUS, docking.GRID_HALF + widest)

    def test_profiles_match_the_untrimmed_receptor(self):
        """The end-to-end guarantee: trimming changes no contact, no count."""
        full = _receptor()
        far = [atom("C", 300 + i, 0, 0, name="C", res_name="LEU", res_seq=50 + i)
               for i in range(20)]
        padded = structure(list(full.protein_atoms) + far)

        lig_atoms = _ligand(3)
        grid = docking.build_grid(padded, (0, 0, 0))
        pose = docking.dock_with_grid(grid, lig_atoms, (0, 0, 0),
                                      seeds=8, mc_steps=5, seed=0)
        component = docking.pose_to_component(pose, "LIG")

        from snaclex import interactions as inter
        before = inter.profile_component(padded, component)
        after = inter.profile_component(
            panel.site_subset(padded, (0, 0, 0)), component
        )
        self.assertEqual(before["counts"], after["counts"])
        self.assertEqual(before["interaction_total"], after["interaction_total"])
        self.assertEqual(before["contact_residue_count"],
                         after["contact_residue_count"])


class TestCofactors(unittest.TestCase):
    """Cofactor-dependent sites must be scorable with the cofactor present."""

    def _structure_with_fad(self):
        from tests.fixtures import component
        prot = [atom("C", 3, 0, 0, name="C", res_name="LEU", res_seq=i)
                for i in range(6)]
        fad = component("FAD", [
            atom("N", 2.0, 1.0, 0.0, hetero=True, res_name="FAD",
                 chain="B", res_seq=500),
            atom("C", 2.0, 2.0, 0.0, hetero=True, res_name="FAD",
                 chain="B", res_seq=500),
        ], chain="B", res_seq=500)
        return structure(prot, [fad])

    def test_declared_cofactor_is_resolved(self):
        s = self._structure_with_fad()
        found, missing = panel.cofactor_atoms(s, ["FAD"])
        self.assertEqual([c.res_name for c in found], ["FAD"])
        self.assertEqual(missing, [])

    def test_absent_cofactor_is_reported_not_raised(self):
        s = self._structure_with_fad()
        found, missing = panel.cofactor_atoms(s, ["FAD", "SAM"])
        self.assertEqual([c.res_name for c in found], ["FAD"])
        self.assertEqual(missing, ["SAM"])

    def test_no_declaration_resolves_to_nothing(self):
        s = self._structure_with_fad()
        self.assertEqual(panel.cofactor_atoms(s, None), ([], []))
        self.assertEqual(panel.cofactor_atoms(s, []), ([], []))

    def test_cofactor_atoms_reach_the_grid_builder(self):
        h = PanelHarness(targets=1)
        h.structures["T0"] = self._structure_with_fad()
        out = h.run(
            [{"id": "T0", "pdb": "T0", "cofactors": ["FAD"]}],
            _ligands(["small"]),
        )
        _pdb, _center, extra = h.grid_builds[0]
        self.assertEqual(len(extra), 2)
        self.assertEqual(out["targets"][0]["cofactors_included"], ["FAD"])
        self.assertNotIn("warning", out["targets"][0])

    def test_missing_cofactor_warns_on_the_target(self):
        h = PanelHarness(targets=1)
        h.structures["T0"] = self._structure_with_fad()
        out = h.run(
            [{"id": "T0", "pdb": "T0", "cofactors": ["FAD", "SAM"]}],
            _ligands(["small"]),
        )
        target = out["targets"][0]
        self.assertEqual(target["cofactors_missing"], ["SAM"])
        self.assertIn("SAM", target["warning"])
        # Still scored — a missing cofactor degrades meaning, not the run.
        self.assertEqual(out["n_scored"], 1)

    def test_cofactor_changes_the_score(self):
        """Proof the cofactor is really in the grid, not just recorded."""
        s = self._structure_with_fad()
        lig = _ligand(2)
        without = docking.dock_with_grid(
            docking.build_grid(s, (0, 0, 0)), lig, (0, 0, 0),
            seeds=12, mc_steps=6, seed=0,
        )
        fad_atoms = [a for c in s.components for a in c.atoms]
        with_cof = docking.dock_with_grid(
            docking.build_grid(s, (0, 0, 0), fad_atoms), lig, (0, 0, 0),
            seeds=12, mc_steps=6, seed=0,
        )
        self.assertNotEqual(without["score"], with_cof["score"])


class TestPanelNormalization(unittest.TestCase):
    def test_z_scores_are_per_target_and_higher_is_better(self):
        h = PanelHarness(targets=2)
        out = h.run(_targets(2), _ligands(["small", "medium", "large"]))
        for target in out["targets"]:
            column = [
                c for c in out["cells"]
                if c["target_id"] == target["id"] and c["error"] is None
            ]
            best = min(column, key=lambda c: c["ligand_efficiency"])
            worst = max(column, key=lambda c: c["ligand_efficiency"])
            # Lower ligand efficiency is a better score => higher z.
            self.assertGreater(best["z_target"], worst["z_target"])
            self.assertEqual(best["rank_in_target"], 1)

    def test_z_scores_are_standardized_within_each_column(self):
        h = PanelHarness(targets=2)
        out = h.run(_targets(2), _ligands(["small", "medium", "large"]))
        for target in out["targets"]:
            zs = [
                c["z_target"] for c in out["cells"]
                if c["target_id"] == target["id"] and c["z_target"] is not None
            ]
            # Mean ~0 by construction; tolerance covers the 3-decimal rounding
            # applied to each published z.
            self.assertAlmostEqual(sum(zs) / len(zs), 0.0, places=2)

    def test_single_ligand_panel_reports_z_unavailable(self):
        """One ligand gives no spread; emit None, not a fake zero."""
        h = PanelHarness(targets=1)
        out = h.run(_targets(1), _ligands(["small"]))
        self.assertFalse(out["normalization"]["per_target"]["T0"]["z_available"])
        self.assertIsNone(out["cells"][0]["z_target"])
        self.assertIsNone(out["selectivity"]["small"]["selectivity_gap"])

    def test_selectivity_picks_the_best_z_and_reports_the_gap(self):
        h = PanelHarness(targets=3)
        out = h.run(_targets(3), _ligands(["small", "medium", "large"]))
        for lid, sel in out["selectivity"].items():
            row = sorted(
                (c for c in out["cells"]
                 if c["ligand_id"] == lid and c["z_target"] is not None),
                key=lambda c: c["z_target"], reverse=True,
            )
            self.assertEqual(sel["best_target_id"], row[0]["target_id"])
            self.assertAlmostEqual(
                sel["selectivity_gap"],
                round(row[0]["z_target"] - row[1]["z_target"], 3),
                places=6,
            )


class TestPanelResilience(unittest.TestCase):
    def test_one_bad_ligand_does_not_abort_the_panel(self):
        h = PanelHarness(targets=2, fail_ligand="medium")
        out = h.run(_targets(2), _ligands(["small", "medium"]))
        self.assertEqual(out["n_cells"], 4)
        self.assertEqual(out["n_failed"], 2)
        bad = [c for c in out["cells"] if c["ligand_id"] == "medium"]
        self.assertTrue(all("ligand unavailable" in c["error"] for c in bad))
        good = [c for c in out["cells"] if c["ligand_id"] == "small"]
        self.assertTrue(all(c["error"] is None for c in good))

    def test_one_bad_target_does_not_abort_the_panel(self):
        h = PanelHarness(targets=2, fail_target="T1")
        out = h.run(_targets(2), _ligands(["small", "medium"]))
        self.assertEqual(out["n_scored"], 2)
        bad = [c for c in out["cells"] if c["target_id"] == "T1"]
        self.assertTrue(all("target unavailable" in c["error"] for c in bad))

    def test_failed_cells_are_excluded_from_normalization(self):
        h = PanelHarness(targets=2, fail_ligand="large")
        out = h.run(_targets(2), _ligands(["small", "medium", "large"]))
        for stats in out["normalization"]["per_target"].values():
            self.assertEqual(stats["n_scored"], 2)


class TestPanelGuards(unittest.TestCase):
    def test_empty_inputs_rejected(self):
        h = PanelHarness(targets=1)
        with self.assertRaises(panel.PanelError):
            h.run([], _ligands(["small"]))
        with self.assertRaises(panel.PanelError):
            h.run(_targets(1), [])

    def test_oversized_panel_rejected_before_any_work(self):
        h = PanelHarness(targets=1)
        many_targets = [{"id": f"T{i}", "pdb": "T0"} for i in range(panel.MAX_TARGETS)]
        many_ligands = [{"id": f"L{i}", "query": "small"}
                        for i in range(panel.MAX_LIGANDS)]
        with self.assertRaises(panel.PanelError) as ctx:
            h.run(many_targets, many_ligands)
        self.assertIn("cells", str(ctx.exception))
        self.assertEqual(h.ligand_loads, [])  # nothing was fetched

    def test_too_many_targets_rejected(self):
        h = PanelHarness(targets=1)
        with self.assertRaises(panel.PanelError):
            h.run([{"id": f"T{i}", "pdb": "T0"} for i in range(panel.MAX_TARGETS + 1)],
                  _ligands(["small"]))


class TestPanelProgressAndMeasured(unittest.TestCase):
    def test_progress_reports_every_cell_in_order(self):
        seen = []
        h = PanelHarness(targets=2)
        h.run(_targets(2), _ligands(["small", "medium"]),
              progress=lambda done, total, label: seen.append((done, total)))
        self.assertEqual(seen, [(1, 4), (2, 4), (3, 4), (4, 4)])

    def test_measured_activity_is_attached_per_cell(self):
        h = PanelHarness(targets=2)
        out = h.run(
            _targets(2), _ligands(["small"]),
            measured_activity=lambda lig, tgt: {"source": "test",
                                                "pair": f"{tgt['id']}/{lig['id']}"},
        )
        for cell in out["cells"]:
            self.assertEqual(
                cell["measured"]["pair"],
                f"{cell['target_id']}/{cell['ligand_id']}",
            )

    def test_measured_lookup_failure_does_not_fail_the_cell(self):
        def boom(lig, tgt):
            raise RuntimeError("chembl down")

        h = PanelHarness(targets=1)
        out = h.run(_targets(1), _ligands(["small", "medium"]),
                    measured_activity=boom)
        self.assertEqual(out["n_failed"], 0)
        self.assertTrue(all(c["measured"] is None for c in out["cells"]))


class TestPanelGeometryAndMethods(unittest.TestCase):
    def test_oversized_ligand_is_flagged_not_silently_clipped(self):
        big = [{"element": "C", "x": 4.0 * i, "y": 0.0, "z": 0.0} for i in range(8)]
        geom = panel._check_box_fit(big)
        self.assertFalse(geom["fits_box"])
        self.assertIn("may not fit", geom["note"])

    def test_normal_ligand_fits(self):
        geom = panel._check_box_fit(_ligand(3))
        self.assertTrue(geom["fits_box"])
        self.assertIsNone(geom["note"])

    def test_methods_block_states_limits_and_search(self):
        h = PanelHarness(targets=2)
        out = h.run(_targets(2), _ligands(["small", "medium"]), seeds=8, mc_steps=5)
        methods = out["methods"]
        self.assertEqual(methods["search"]["seeds"], 8)
        self.assertEqual(methods["n_targets"], 2)
        self.assertTrue(methods["limitations"])
        self.assertIn("not an affinity", methods["scoring"].lower())

    def test_normalization_note_warns_raw_scores_are_not_comparable(self):
        h = PanelHarness(targets=2)
        out = h.run(_targets(2), _ligands(["small", "medium"]))
        self.assertIn("NOT comparable", out["normalization"]["note"])


if __name__ == "__main__":
    unittest.main()
