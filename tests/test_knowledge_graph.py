"""Tests for the provenance-aware knowledge-graph projection. Pure, offline."""

import unittest

from snaclex import evidence, knowledge_graph


def _record():
    return {
        "canonical_id": "SNX:PRT:9606:P00519:1",
        "preferred_name": "Tyrosine-protein kinase ABL1",
        "gene": {"symbol": "ABL1", "synonyms": []},
        "organism": {"scientific_name": "Homo sapiens", "taxon_id": 9606},
        "review_status": "reviewed",
        "accessions": {"uniprot_primary": "P00519", "uniparc": "UPI000012ABCD"},
        "sequence": {"length": 1130, "checksum": {"crc64": "A1B2C3D4E5F60718"}},
        "isoforms": [{"id": "P00519-2", "name": "IB"}],
        "domains_motifs": [{"name": "Protein kinase", "start": 242, "end": 493}],
        "variants": [{"position": 315, "wt": "T", "mut": "I", "description": "CML"}],
        "disease_associations": [{"name": "Leukemia", "accession": "DI-001"}],
        "literature": [{"pmid": "3018722", "title": "Cloning", "source": "UniProtKB reference"}],
        "structures": {
            "experimental": [{"pdb_id": "1IEP"}],
            "homologous": [{"pdb_id": "2GQG", "sequence_identity": 0.62}],
            "predicted": [{"model_id": "AF-P00519-F1", "global_plddt_mean": 88.0}],
        },
        "xref_graph": {
            "nodes": [{"ns": "UniProt", "id": "P00519"},
                      {"ns": "RefSeq", "id": "NP_005148.2"}],
            "edges": [{"from": "UniProt:P00519", "to": "RefSeq:NP_005148.2",
                       "rel": "maps_to", "source": "UniProtKB"}],
        },
        "provenance_summary": {"retrieved_utc": "2026-07-12T00:00:00Z"},
    }


class TestBuildGraph(unittest.TestCase):
    def setUp(self):
        self.g = knowledge_graph.build_graph(_record())

    def test_core_nodes_present(self):
        types = {n["type"] for n in self.g["nodes"]}
        for t in ("Protein", "Gene", "Organism", "ProteinSequence", "Isoform",
                  "Domain", "Variant", "Disease", "Publication", "Structure",
                  "DatabaseRecord"):
            self.assertIn(t, types)

    def test_structure_edges_typed(self):
        rels = {e["rel"] for e in self.g["edges"]}
        self.assertIn("has_structure", rels)
        self.assertIn("has_predicted_structure", rels)
        # Predicted edge is flagged predicted, not experimental.
        pred = [e for e in self.g["edges"] if e["rel"] == "has_predicted_structure"][0]
        self.assertEqual(pred["provenance"]["experimental_or_predicted"], "predicted")
        self.assertEqual(pred["provenance"]["source_database"], "AlphaFold DB")

    def test_homolog_structure_flagged(self):
        homolog = [e for e in self.g["edges"]
                   if e["to"] == "PDB:2GQG"][0]
        self.assertEqual(homolog["provenance"]["evidence_type"], "homology-transferred")

    def test_every_edge_has_provenance(self):
        for e in self.g["edges"]:
            p = e["provenance"]
            self.assertIn("source_database", p)
            self.assertIn("evidence_type", p)
            self.assertIn("experimental_or_predicted", p)
            self.assertIn("retrieval_date", p)

    def test_encoded_by_edge(self):
        e = [x for x in self.g["edges"] if x["rel"] == "encoded_by"][0]
        self.assertEqual(e["to"], "GENE:ABL1")

    def test_counts(self):
        self.assertEqual(self.g["counts"]["nodes"], len(self.g["nodes"]))
        self.assertEqual(self.g["counts"]["edges"], len(self.g["edges"]))
        self.assertIn("has_structure", self.g["counts"]["by_edge_type"])


class TestEvidenceInGraph(unittest.TestCase):
    def test_typed_evidence_edge_and_software(self):
        ev = [
            evidence.docking_evidence(
                protein_ref="SNX:PRT:9606:P00519:1",
                chemical_ref="PubChem:CID:5291", score=-7.2,
                software="SnaCleX docker", params={}, seed=42),
        ]
        rec = _record()
        g = knowledge_graph.build_graph(rec, evidence=ev)
        binds = [e for e in g["edges"] if e["rel"] == "binds"]
        self.assertTrue(binds)
        self.assertIn("level E", binds[0]["provenance"]["evidence_type"])
        self.assertEqual(binds[0]["provenance"]["experimental_or_predicted"],
                         "predicted/inferred")
        self.assertEqual(binds[0]["provenance"]["software_version"], "SnaCleX docker")
        self.assertIn("Compound", {n["type"] for n in g["nodes"]})


if __name__ == "__main__":
    unittest.main()
