"""Tests for research-domain workspace lenses (Phase 3). Pure, offline."""

import unittest

from snaclex import workspaces


def _rec(gene=None, name="", keywords=(), domains=(), variants=(),
         diseases=(), isoforms=(), synonyms=()):
    return {
        "gene": {"symbol": gene, "synonyms": list(synonyms)},
        "preferred_name": name,
        "keywords": [{"name": k, "category": "x"} for k in keywords],
        "domains_motifs": [{"name": d} for d in domains],
        "variants": list(variants),
        "disease_associations": list(diseases),
        "isoforms": list(isoforms),
        "accessions": {"uniprot_primary": "P00000"},
    }


class TestImmunology(unittest.TestCase):
    def test_cytokine_by_keyword(self):
        v = workspaces.immunology_view(_rec(gene="IL6", name="Interleukin-6",
                                            keywords=["Cytokine"]))
        self.assertTrue(v["in_scope"])
        self.assertIn("cytokine / cytokine receptor", v["categories"])

    def test_chemokine_by_gene(self):
        v = workspaces.immunology_view(_rec(gene="CXCL8", name="Interleukin-8"))
        self.assertIn("chemokine", v["categories"])

    def test_checkpoint(self):
        v = workspaces.immunology_view(_rec(gene="PDCD1", name="Programmed cell death protein 1"))
        self.assertIn("immune checkpoint", v["categories"])

    def test_hla_not_collapsed_and_allele_preserved(self):
        v = workspaces.immunology_view(
            _rec(gene="HLA-A", name="HLA class I histocompatibility antigen",
                 keywords=["MHC I"]),
            query="HLA-A*02:01")
        self.assertTrue(v["hla_mhc"]["is_hla_mhc"])
        self.assertEqual(v["hla_mhc"]["allele_in_query"], "HLA-A*02:01")
        self.assertIn("NOT collapsed", v["hla_mhc"]["note"])

    def test_out_of_scope(self):
        v = workspaces.immunology_view(_rec(gene="ALB", name="Albumin"))
        self.assertFalse(v["in_scope"])


class TestOncology(unittest.TestCase):
    def test_oncogene_by_keyword(self):
        v = workspaces.oncology_view(_rec(gene="ABL1", keywords=["Proto-oncogene"]))
        roles = {r["role"] for r in v["roles"]}
        self.assertIn("oncogene", roles)
        self.assertTrue(any("curated" in r["basis"] for r in v["roles"]))
        self.assertTrue(v["cancer_associated"])

    def test_tumor_suppressor_and_dna_repair(self):
        v = workspaces.oncology_view(
            _rec(gene="BRCA1", keywords=["Tumor suppressor", "DNA repair"]))
        roles = {r["role"] for r in v["roles"]}
        self.assertIn("tumor suppressor", roles)
        self.assertIn("DNA repair", roles)

    def test_kinase_by_domain(self):
        v = workspaces.oncology_view(_rec(gene="XYZ", domains=["Protein kinase"]))
        self.assertIn("kinase", {r["role"] for r in v["roles"]})

    def test_cancer_variants_and_separation_note(self):
        v = workspaces.oncology_view(_rec(
            gene="BRAF", keywords=["Proto-oncogene"],
            variants=[{"position": 600, "wt": "V", "mut": "E", "description": "in melanoma"}],
            diseases=[{"name": "Colorectal cancer", "description": "somatic"}]))
        self.assertTrue(v["cancer_associated_variants"])
        self.assertTrue(v["cancer_disease_associations"])
        self.assertIn("NO claim", v["separation_note"])

    def test_seed_list_flagged_heuristic(self):
        v = workspaces.oncology_view(_rec(gene="KRAS"))  # no keyword, seed only
        self.assertTrue(any("heuristic" in r["basis"] for r in v["roles"]))


class TestGenetics(unittest.TestCase):
    def test_counts_and_clinical_note(self):
        v = workspaces.genetics_view(_rec(
            gene="TP53", synonyms=["P53"],
            isoforms=[{"id": "P04637-1"}, {"id": "P04637-2"}],
            variants=[{"position": 175, "description": "in cancer"}]))
        self.assertEqual(v["isoform_count"], 2)
        self.assertEqual(v["variant_count"], 1)
        self.assertIn("never as a clinical recommendation", v["clinical_note"])


class TestApplyViews(unittest.TestCase):
    def test_selective(self):
        out = workspaces.apply_views(_rec(gene="ABL1"), {"oncology"})
        self.assertIn("oncology", out)
        self.assertNotIn("immunology", out)


if __name__ == "__main__":
    unittest.main()
