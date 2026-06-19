"""Tests for the method-transparency / provenance blocks (Phase 2)."""

import json
import unittest

from snaclex import __version__, evolution, pockets, provenance


class TestProvenance(unittest.TestCase):
    def test_tool_carries_version(self):
        self.assertIn(__version__, provenance.tool())

    def test_pocket_methods_shape(self):
        m = provenance.pocket_methods()
        for key in ("tool", "method", "method_family", "parameters",
                    "scoring", "interpretation", "limitations", "disclaimer"):
            self.assertIn(key, m)
        self.assertTrue(m["limitations"])
        self.assertIn("Research-only", m["disclaimer"])

    def test_pocket_methods_pull_real_constants(self):
        m = provenance.pocket_methods()
        self.assertEqual(m["parameters"]["psp_threshold"], pockets.PSP_THRESHOLD)
        self.assertEqual(m["parameters"]["scan_range_steps"], pockets.SCAN_RANGE)
        self.assertEqual(m["parameters"]["min_pocket_points"], pockets.MIN_POCKET_POINTS)

    def test_evolution_methods_shape(self):
        m = provenance.evolution_methods()
        for key in ("tool", "method", "parameters", "interpretation",
                    "limitations", "disclaimer"):
            self.assertIn(key, m)
        self.assertEqual(
            m["parameters"]["max_alignment_sequences"], evolution.MAX_ALIGN_SEQS
        )

    def test_new_method_blocks_shape(self):
        for fn in (provenance.interface_methods, provenance.variant_methods,
                   provenance.hla_methods, provenance.esm_methods):
            m = fn()
            for key in ("tool", "method", "method_family", "parameters",
                        "interpretation", "limitations", "disclaimer"):
                self.assertIn(key, m)
            self.assertTrue(m["limitations"])
            self.assertIn("Research-only", m["disclaimer"])

    def test_blocks_are_json_serializable(self):
        # They are embedded directly in API JSON responses.
        json.dumps(provenance.pocket_methods())
        json.dumps(provenance.evolution_methods())
        json.dumps(provenance.interface_methods())
        json.dumps(provenance.variant_methods())
        json.dumps(provenance.hla_methods())
        json.dumps(provenance.esm_methods())


if __name__ == "__main__":
    unittest.main()
