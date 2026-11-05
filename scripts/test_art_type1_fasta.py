import unittest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from audit_art_type1_fasta import six_frame_orfs, protein_interval_to_source, matching_annotations


class TestType1RescueCoordinates(unittest.TestCase):
    def test_six_frame_translation_and_reverse_mapping(self):
        dna = Seq("GGTAAATG" + "GCT" * 45 + "TAACCG" + "NNN" + "CAT" * 42 + "TGA")
        for source in (dna, dna.reverse_complement()):
            orfs = list(six_frame_orfs(source))
            self.assertTrue(orfs)
            for row in orfs:
                fragment = source[row["start0"]:row["end0"]]
                if row["strand"] == -1:
                    fragment = fragment.reverse_complement()
                self.assertEqual(str(fragment.translate(table=11)), row["peptide"])
                row["length_aa"] = len(row["peptide"])
                a, b = protein_interval_to_source(row, 2, 12)
                mapped = source[a:b]
                if row["strand"] == -1:
                    mapped = mapped.reverse_complement()
                self.assertEqual(str(mapped.translate(table=11)), row["peptide"][2:12])

    def test_ambiguity_and_stops_do_not_join_fragments(self):
        dna = "GCT" * 41 + "NNN" + "GCT" * 43 + "TAA" + "GCT" * 44
        rows = [r for r in six_frame_orfs(dna) if r["strand"] == 1
                and r["frame_in_oriented_sequence"] == 0]
        self.assertEqual([len(r["peptide"]) for r in rows], [41, 43, 44])
        self.assertEqual(rows[0]["boundary_3"], "ambiguous")
        self.assertEqual(rows[1]["boundary_3"], "stop")

    def test_partial_non_start_fragments_are_retained(self):
        rows = list(six_frame_orfs("GCT" * 42))
        self.assertIn("A" * 42, [r["peptide"] for r in rows])

    def test_annotation_overlap_requires_strand_and_frame(self):
        f = SeqFeature(SimpleLocation(101, 404, strand=-1), type="CDS",
                       qualifiers={"protein_id": ["control"]})
        r = dict(start0=101, end0=434, strand=-1)
        self.assertEqual(matching_annotations(r, [f]), "control")
        self.assertEqual(matching_annotations(dict(r, end0=435), [f]), "")
        self.assertEqual(matching_annotations(dict(r, strand=1), [f]), "")


if __name__ == "__main__":
    unittest.main()
