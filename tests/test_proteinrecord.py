"""Tests for the unified ProteinRecord assembly (pure, offline)."""

import unittest

from snaclex import proteinrecord, uniprot

_ENTRY = {
    "primaryAccession": "P00519", "uniProtkbId": "ABL1_HUMAN",
    "entryType": "UniProtKB reviewed (Swiss-Prot)", "entryAudit": {"entryVersion": 1},
    "sequence": {"value": "MLEICLKLVGCKSKKGLSSS", "length": 20, "crc64": "AABB", "md5": "cc"},
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "genes": [{"geneName": {"value": "ABL1"}}],
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Tyrosine-protein kinase ABL1"}}},
    "comments": [], "features": [
        {"type": "Domain", "location": {"start": {"value": 1}, "end": {"value": 10}},
         "description": "Protein kinase"}],
    "uniProtKBCrossReferences": [{"database": "PDB", "id": "1IEP", "properties": []}],
}


def _resolution(**over):
    base = {
        "query": "ABL1", "query_type": "uniprot_accession",
        "chosen": {"accession": "P00519", "taxon_id": 9606, "review_status": "reviewed"},
        "xref_graph": {"nodes": [], "edges": []},
        "identity_key": {"crc64": "AABB"}, "warnings": [],
        "_uniprot_fragment": uniprot.parse_entry(_ENTRY)["data"],
        "_uniparc_fragment": None,
    }
    base.update(over)
    return base


class TestBuildRecord(unittest.TestCase):
    def test_identity_and_fields(self):
        rec = proteinrecord.build_record(_resolution())
        self.assertEqual(rec["canonical_id"], "SNX:PRT:9606:P00519:1")
        self.assertEqual(rec["review_status"], "reviewed")
        self.assertEqual(rec["preferred_name"], "Tyrosine-protein kinase ABL1")
        self.assertEqual(rec["accessions"]["uniprot_primary"], "P00519")
        self.assertTrue(rec["domains_motifs"])
        self.assertIn("UniProtKB", rec["provenance_summary"]["sources_used"])

    def test_sequence_analysis_always_runs(self):
        rec = proteinrecord.build_record(_resolution())
        self.assertIsNotNone(rec["sequence_analysis"])
        self.assertEqual(rec["sequence_analysis"]["length"], 20)
        self.assertIn("isoelectric_point", rec["sequence_analysis"])

    def test_checksum_filled_when_absent(self):
        entry = dict(_ENTRY)
        entry["sequence"] = {"value": "MLEICLKLVGCKSKKGLSSS", "length": 20}
        res = _resolution(_uniprot_fragment=uniprot.parse_entry(entry)["data"])
        rec = proteinrecord.build_record(res)
        self.assertIn("sha256", rec["sequence"]["checksum"])

    def test_sequence_only_record(self):
        res = {
            "query": ">x\nMLEICLKLVGCKSKKGLSSS", "query_type": "raw_sequence",
            "sequence": "MLEICLKLVGCKSKKGLSSS", "identity_key": {"crc64": "AABB"},
            "warnings": ["novel/unmatched sequence — no UniProt/UniParc record; sequence-only"],
            "xref_graph": {"nodes": [], "edges": []},
        }
        rec = proteinrecord.build_record(res)
        self.assertEqual(rec["review_status"], "sequence-only")
        self.assertIsNone(rec["accessions"]["uniprot_primary"])
        self.assertEqual(rec["canonical_id"], "SNX:PRT:NA:SEQONLY:1")
        self.assertIsNotNone(rec["sequence_analysis"])
        self.assertTrue(any("sequence-only" in w for w in rec["warnings"]))

    def test_unreviewed_warning(self):
        entry = dict(_ENTRY)
        entry["entryType"] = "UniProtKB unreviewed (TrEMBL)"
        res = _resolution(_uniprot_fragment=uniprot.parse_entry(entry)["data"])
        rec = proteinrecord.build_record(res)
        self.assertEqual(rec["review_status"], "unreviewed")
        self.assertTrue(any("unreviewed" in w for w in rec["warnings"]))


class TestAttachStructures(unittest.TestCase):
    def test_predicted_tier_flags_warning(self):
        rec = proteinrecord.build_record(_resolution())
        availability = {"experimental": [], "homologous": [],
                        "predicted": [{"model_id": "AF-P00519-F1", "docking_suitable": False}],
                        "tier": "predicted", "best_for_docking": None, "warnings": []}
        proteinrecord.attach_structures(rec, availability)
        self.assertEqual(rec["structures"]["tier"], "predicted")
        self.assertTrue(any("predicted model" in w for w in rec["warnings"]))


if __name__ == "__main__":
    unittest.main()
