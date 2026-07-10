"""Offline parse tests for the new source adapters.

Each adapter keeps its network entry point thin and its parsing pure; these
tests exercise the parsing against realistically shaped fixture payloads, so no
live UniProt/UniParc/NCBI/AlphaFold/RCSB/PubChem call is made. Fixtures mirror
the documented response shapes (see docs/platform/03-source-integration.md).
"""

import unittest

from snaclex import alphafold, ncbi, pubchem, rcsb, uniparc, uniprot


# --- UniProtKB --------------------------------------------------------------
UNIPROT_ENTRY = {
    "primaryAccession": "P00519",
    "uniProtkbId": "ABL1_HUMAN",
    "entryType": "UniProtKB reviewed (Swiss-Prot)",
    "entryAudit": {"entryVersion": 275},
    "sequence": {"value": "MLEICLKLVG", "length": 1130, "molWeight": 122873,
                 "crc64": "A1B2C3D4E5F60718", "md5": "deadbeef"},
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "genes": [{"geneName": {"value": "ABL1"},
               "synonyms": [{"value": "ABL"}, {"value": "JTK7"}]}],
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Tyrosine-protein kinase ABL1"}}},
    "comments": [
        {"commentType": "FUNCTION", "texts": [{"value": "Non-receptor tyrosine-protein kinase."}]},
        {"commentType": "SUBCELLULAR LOCATION",
         "subcellularLocations": [{"location": {"value": "Cytoplasm"}}]},
        {"commentType": "ALTERNATIVE PRODUCTS",
         "isoforms": [{"isoformIds": ["P00519-1"], "name": {"value": "IB"},
                       "isoformSequenceStatus": "Displayed"}]},
    ],
    "features": [
        {"type": "Domain", "location": {"start": {"value": 242}, "end": {"value": 493}},
         "description": "Protein kinase"},
        {"type": "Active site", "location": {"start": {"value": 363}, "end": {"value": 363}},
         "description": "Proton acceptor"},
        {"type": "Natural variant", "location": {"start": {"value": 315}, "end": {"value": 315}},
         "description": "in CML", "alternativeSequence": {"originalSequence": "T",
                                                          "alternativeSequences": ["I"]}},
        {"type": "Modified residue", "location": {"start": {"value": 393}, "end": {"value": 393}},
         "description": "Phosphotyrosine"},
    ],
    "uniProtKBCrossReferences": [
        {"database": "RefSeq", "id": "NP_005148.2", "properties": []},
        {"database": "PDB", "id": "1IEP", "properties": []},
        {"database": "GeneID", "id": "25", "properties": []},
        {"database": "AlphaFoldDB", "id": "P00519", "properties": []},
    ],
}


class TestUniprot(unittest.TestCase):
    def test_parse_entry(self):
        frag = uniprot.parse_entry(UNIPROT_ENTRY)
        d = frag["data"]
        self.assertEqual(d["accession"], "P00519")
        self.assertEqual(d["review_status"], "reviewed")
        self.assertEqual(d["gene"]["symbol"], "ABL1")
        self.assertIn("ABL", d["gene"]["synonyms"])
        self.assertEqual(d["preferred_name"], "Tyrosine-protein kinase ABL1")
        self.assertEqual(d["sequence"]["length"], 1130)
        self.assertEqual(frag["provenance"]["source"], "UniProtKB")

    def test_features_classified(self):
        d = uniprot.parse_entry(UNIPROT_ENTRY)["data"]
        self.assertTrue(any(x["name"] == "Protein kinase" for x in d["domains_motifs"]))
        self.assertTrue(any(x["kind"] == "active_site" for x in d["residue_annotations"]))
        var = d["variants"][0]
        self.assertEqual((var["wt"], var["position"], var["mut"]), ("T", 315, "I"))
        self.assertTrue(d["ptms"])

    def test_cross_references(self):
        d = uniprot.parse_entry(UNIPROT_ENTRY)["data"]
        self.assertIn("RefSeq", d["cross_references"])
        self.assertIn("GeneID", d["cross_references"])

    def test_reviewed_vs_unreviewed(self):
        self.assertEqual(uniprot.review_status("UniProtKB reviewed (Swiss-Prot)"), "reviewed")
        self.assertEqual(uniprot.review_status("UniProtKB unreviewed (TrEMBL)"), "unreviewed")

    def test_parse_search(self):
        data = {"results": [{
            "primaryAccession": "P00519", "uniProtkbId": "ABL1_HUMAN",
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "proteinDescription": {"recommendedName": {"fullName": {"value": "ABL1"}}},
            "genes": [{"geneName": {"value": "ABL1"}}],
            "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
            "sequence": {"length": 1130}}]}
        cands = uniprot.parse_search(data)
        self.assertEqual(cands[0]["accession"], "P00519")
        self.assertEqual(cands[0]["review_status"], "reviewed")


# --- UniParc ----------------------------------------------------------------
class TestUniparc(unittest.TestCase):
    ENTRY = {
        "uniParcId": "UPI00001AB2CD",
        "sequence": {"value": "MLEICLKLVG", "length": 1130, "crc64": "A1B2", "md5": "x"},
        "uniParcCrossReferences": [
            {"database": "UniProtKB/Swiss-Prot", "id": "P00519", "active": True,
             "version": 4, "versionI": 1, "created": "2000-01-01", "lastUpdated": "2020-01-01"},
            {"database": "RefSeq", "id": "NP_005148.2", "active": True},
            {"database": "EMBL", "id": "OLD123", "active": False},
        ],
    }

    def test_parse_entry(self):
        frag = uniparc.parse_entry(self.ENTRY)
        d = frag["data"]
        self.assertEqual(d["upi"], "UPI00001AB2CD")
        self.assertIn("UniProtKB/Swiss-Prot", d["active_databases"])
        self.assertEqual(len(d["obsolete_references"]), 1)

    def test_uniprot_accessions(self):
        d = uniparc.parse_entry(self.ENTRY)["data"]
        self.assertEqual(uniparc.uniprot_accessions(d), ["P00519"])


# --- NCBI -------------------------------------------------------------------
class TestNcbi(unittest.TestCase):
    def test_parse_fasta(self):
        fasta = (">NP_005148.2 tyrosine-protein kinase ABL1 isoform a [Homo sapiens]\n"
                 "MLEICLKLVG\nCKSKKGLSSS\n")
        rec = ncbi.parse_fasta(fasta)
        self.assertEqual(rec["accession"], "NP_005148.2")
        self.assertEqual(rec["organism"], "Homo sapiens")
        self.assertEqual(rec["sequence"], "MLEICLKLVGCKSKKGLSSS")

    def test_parse_esummary(self):
        data = {"result": {"uids": ["25"], "25": {
            "title": "ABL1", "organism": "Homo sapiens", "taxid": 9606,
            "slen": 1130, "accessionversion": "NP_005148.2", "gi": "1"}}}
        s = ncbi.parse_esummary(data)
        self.assertEqual(s["taxid"], 9606)
        self.assertEqual(s["length"], 1130)

    def test_refseq_class(self):
        self.assertEqual(ncbi.refseq_class("NP_005148.2"), "refseq-curated")
        self.assertEqual(ncbi.refseq_class("XP_011529125"), "refseq-predicted")
        self.assertEqual(ncbi.refseq_class("WP_000000001"), "refseq-nonredundant")
        self.assertIsNone(ncbi.refseq_class("P00519"))

    def test_build_fragment_predicted_flag(self):
        frag = ncbi.build_fragment(
            "XP_011529125", {"accession": "XP_011529125", "sequence": "MLE"}, {})
        self.assertEqual(frag["provenance"]["evidence_category"], "predicted")


# --- AlphaFold --------------------------------------------------------------
def _ca_line(serial, resseq, b):
    line = list(" " * 80)

    def put(col, s):
        for i, ch in enumerate(s):
            line[col - 1 + i] = ch

    put(1, "ATOM")
    put(7, f"{serial:>5}")
    put(13, "CA")
    put(18, "ALA")
    put(22, "A")
    put(23, f"{resseq:>4}")
    put(31, f"{0.0:>8.3f}")
    put(39, f"{0.0:>8.3f}")
    put(47, f"{0.0:>8.3f}")
    put(55, f"{1.0:>6.2f}")
    put(61, f"{b:>6.2f}")
    put(77, " C")
    return "".join(line)


class TestAlphafold(unittest.TestCase):
    def test_parse_prediction(self):
        data = [{"entryId": "AF-P00519-F1", "uniprotAccession": "P00519",
                 "uniprotStart": 1, "uniprotEnd": 10, "uniprotSequence": "MLEICLKLVG",
                 "latestVersion": 4, "modelCreatedDate": "2022-06-01",
                 "globalMetricValue": 71.4,
                 "pdbUrl": "https://x/AF-P00519-F1-model_v4.pdb",
                 "paeImageUrl": "https://x/pae.png"}]
        models = alphafold.parse_prediction(data)
        m = models[0]["data"]
        self.assertEqual(m["model_id"], "AF-P00519-F1")
        self.assertEqual(m["coverage_fraction"], 1.0)
        self.assertTrue(m["pae_available"])
        self.assertEqual(models[0]["provenance"]["evidence_category"], "predicted")

    def test_parse_plddt_from_bfactor(self):
        pdb = "\n".join([_ca_line(1, 1, 95.0), _ca_line(2, 2, 40.0)])
        plddt = alphafold.parse_plddt(pdb)
        self.assertEqual(plddt, [95.0, 40.0])

    def test_confidence_bands(self):
        bands = alphafold.confidence_bands([95.0, 40.0])
        self.assertEqual(bands["very_high"], 0.5)
        self.assertEqual(bands["very_low"], 0.5)

    def test_low_confidence_regions(self):
        regions = alphafold.low_confidence_regions([95, 40, 40, 95], threshold=70, min_len=1)
        self.assertEqual(regions, [{"start": 2, "end": 3}])

    def test_docking_gate(self):
        self.assertTrue(alphafold.docking_assessment({"global_plddt_mean": 95})["docking_suitable"])
        self.assertFalse(alphafold.docking_assessment({"global_plddt_mean": 40})["docking_suitable"])

    def test_build_confidence_high(self):
        pdb = "\n".join(_ca_line(i, i, 95.0) for i in range(1, 6))
        conf = alphafold.build_confidence({"model_id": "AF-X-F1"}, pdb)
        self.assertEqual(conf["n_residues"], 5)
        self.assertEqual(conf["mean_plddt"], 95.0)
        self.assertTrue(conf["docking_suitable"])
        self.assertEqual(conf["low_confidence_regions"], [])

    def test_build_confidence_low_blocks_docking(self):
        pdb = "\n".join(_ca_line(i, i, 40.0) for i in range(1, 6))
        conf = alphafold.build_confidence({"model_id": "AF-X-F1"}, pdb)
        self.assertFalse(conf["docking_suitable"])
        self.assertTrue(conf["low_confidence_regions"])
        self.assertTrue(conf["disordered_regions"])


# --- RCSB search parse ------------------------------------------------------
class TestRcsbSearch(unittest.TestCase):
    def test_parse_polymer_entity_search(self):
        data = {"result_set": [{"identifier": "1IEP_1"}, {"identifier": "2HYY_1"}]}
        out = rcsb.parse_polymer_entity_search(data)
        self.assertEqual(out[0]["pdb_id"], "1IEP")
        self.assertEqual(out[0]["entity"], "1")
        self.assertEqual({o["pdb_id"] for o in out}, {"1IEP", "2HYY"})

    def test_parse_sequence_search(self):
        data = {"result_set": [
            {"identifier": "2GQG_1", "services": [{"nodes": [{"match_context": [
                {"sequence_identity": 0.62, "evalue": 1e-40}]}]}]},
            {"identifier": "3XYZ_1", "services": [{"nodes": [{"match_context": [
                {"sequence_identity": 0.25, "evalue": 0.1}]}]}]},
        ]}
        out = rcsb.parse_sequence_search(data)
        self.assertEqual(out[0]["pdb_id"], "2GQG")
        self.assertEqual(out[0]["identity"], 0.62)
        self.assertAlmostEqual(out[0]["evalue"], 1e-40)

    def test_sequence_search_rejects_short(self):
        # No network call for a too-short sequence.
        self.assertEqual(rcsb.sequence_search("ACDEFG"), [])


# --- PubChem BioAssay parse -------------------------------------------------
class TestPubchemAssay(unittest.TestCase):
    def test_parse_assay_summary(self):
        data = {"Table": {
            "Columns": {"Column": ["AID", "CID", "Activity Outcome", "Target GeneID",
                                   "Activity Name", "Activity Value [uM]", "Assay Name",
                                   "Assay Type", "PubMed ID"]},
            "Row": [
                {"Cell": ["372", "5291", "Active", "25", "IC50", "0.038",
                          "ABL1 assay", "confirmatory", "12345"]},
                {"Cell": ["999", "5291", "Inactive", "25", "IC50", "",
                          "other", "screening", ""]},
            ]}}
        rows = pubchem.parse_assay_summary(data)
        self.assertEqual(rows[0]["aid"], "372")
        self.assertEqual(rows[0]["outcome"], "Active")
        self.assertEqual(rows[0]["target_geneid"], "25")
        self.assertAlmostEqual(rows[0]["activity_value_uM"], 0.038)
        self.assertEqual(rows[0]["pmid"], "12345")
        self.assertEqual(rows[1]["outcome"], "Inactive")  # negatives kept
        self.assertIsNone(rows[1]["activity_value_uM"])


if __name__ == "__main__":
    unittest.main()
