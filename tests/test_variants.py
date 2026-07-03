"""Offline tests for the genome variant bridge (parsing + structure mapping)."""

import unittest

from snaclex import variants
from snaclex.pdbparse import Atom
from tests.fixtures import structure as make_structure

_ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
    "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
    "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
    "Y": "TYR", "V": "VAL",
}

# A short synthetic "protein" sequence used as both UniProt and structure seq.
SEQ = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"


def _structure_from_seq(seq, chain="A", start=1):
    atoms = [
        Atom(i + 1, "CA", _ONE_TO_THREE[a], chain, start + i, "", i * 3.8, 0, 0, "C", False)
        for i, a in enumerate(seq)
    ]
    return make_structure(protein=atoms)


class TestVariantParsing(unittest.TestCase):
    def test_parse_missense_only(self):
        features = [
            {  # missense with ClinVar + gnomAD
                "type": "VARIANT", "begin": "10", "end": "10",
                "wildType": "A", "alternativeSequence": "T",
                "clinicalSignificances": [{"type": "Pathogenic"}],
                "populationFrequencies": [{"source": "gnomAD", "frequency": 0.0004}],
                "xrefs": [{"name": "dbSNP", "id": "rs1"}, {"name": "ClinVar", "id": "RCV1"}],
            },
            {  # deletion (multi-residue) — must be skipped
                "type": "VARIANT", "begin": "20", "end": "22",
                "wildType": "AAA", "alternativeSequence": "-",
            },
            {  # no ClinVar and no gnomAD — skipped as noise
                "type": "VARIANT", "begin": "30", "end": "30",
                "wildType": "K", "alternativeSequence": "R",
            },
        ]
        parsed = variants._parse_variants(features)
        self.assertEqual(len(parsed), 1)
        v = parsed[0]
        self.assertEqual((v["ref"], v["alt"], v["position"]), ("A", "T", 10))
        self.assertEqual(v["clinical_significance"], "pathogenic")
        self.assertAlmostEqual(v["allele_frequency"], 0.0004)
        self.assertEqual(v["clinvar"], "RCV1")

    def test_gnomad_af_takes_max(self):
        af = variants._gnomad_af([
            {"source": "gnomAD", "frequency": 0.001},
            {"source": "gnomAD exomes", "frequency": 0.01},
            {"source": "1000Genomes", "frequency": 0.9},
        ])
        self.assertAlmostEqual(af, 0.01)

    def test_canonical_sig_picks_worst(self):
        rank, label = variants._canonical_sig(
            [{"type": "Benign"}, {"type": "Likely pathogenic"}]
        )
        self.assertEqual(label, "likely pathogenic")


class TestVariantMapping(unittest.TestCase):
    def test_map_positions_identity(self):
        struct = _structure_from_seq(SEQ, chain="A", start=1)
        pos_map, mapped = variants._map_positions(struct, SEQ)
        # UniProt position 10 should map to structure residue res_seq 10 on A.
        self.assertIn(10, pos_map)
        self.assertEqual(pos_map[10][0][0], "A")
        self.assertEqual(pos_map[10][0][1], 10)
        self.assertGreater(mapped, 0)

    def test_map_positions_offset_numbering(self):
        # Structure numbered from 100 (author numbering) still maps by sequence.
        struct = _structure_from_seq(SEQ, chain="B", start=100)
        pos_map, _ = variants._map_positions(struct, SEQ)
        self.assertEqual(pos_map[5][0], ("B", 104, _ONE_TO_THREE[SEQ[4]]))

    def test_analyze_end_to_end_with_monkeypatched_fetch(self):
        struct = _structure_from_seq(SEQ, chain="A", start=1)
        orig = variants.fetch_variation
        variants.fetch_variation = lambda acc: {
            "accession": acc, "gene": "TESTG", "sequence": SEQ,
            "variants": [
                {"position": 10, "ref": "A", "alt": "T",
                 "clinical_significance": "pathogenic", "sig_rank": 5,
                 "allele_frequency": 0.0004, "dbsnp": "rs1", "clinvar": "RCV1"},
                {"position": 5, "ref": SEQ[4], "alt": "W",
                 "clinical_significance": "benign", "sig_rank": 0,
                 "allele_frequency": 0.2, "dbsnp": None, "clinvar": None},
            ],
        }
        orig_locus = variants.locus_info
        variants.locus_info = lambda gene: {"gene": gene, "chromosome": "1"}
        try:
            result = variants.analyze(struct, ["P00000"])
        finally:
            variants.fetch_variation = orig
            variants.locus_info = orig_locus
        self.assertTrue(result["available"])
        self.assertEqual(result["uniprot"], "P00000")
        self.assertEqual(result["mapped_variant_count"], 2)
        # Pathogenic residue ranks first.
        self.assertEqual(result["residues"][0]["pathogenicity"], "pathogenic")
        self.assertEqual(result["residues"][0]["res_seq"], 10)

    def test_analyze_unavailable_without_uniprot(self):
        struct = _structure_from_seq(SEQ)
        result = variants.analyze(struct, [])
        self.assertFalse(result["available"])
        self.assertIn("UniProt", result["reason"])


if __name__ == "__main__":
    unittest.main()
