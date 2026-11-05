"""Regression tests for benchmark-exposed alignment failures and locus evidence gates."""
import json
from pathlib import Path
import tempfile
import unittest

from Bio import SeqIO
from art_array_detector import ROOT, align_boundaries
from art_locus_evidence import classify_locus


class ValidationTests(unittest.TestCase):
    def test_chance_first_seed_does_not_veto_majority_repeat_alignment(self):
        # Frozen independently planted motif, with one chance approximate seed
        # before six true copies. The old first-copy anchor returned no block.
        record=SeqIO.read(ROOT/'results/art_validation/loci/long_spacer_120_1.gb','genbank')
        cfg=json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp') as temp:
            out=Path(temp);(out/'alignments').mkdir()
            result=align_boundaries(str(record.seq[:6000]),[3798,4006,4154,4302,4450,4598,4746],10,cfg,out,'outlier')
        self.assertEqual(result['status'],'aligned')
        self.assertEqual(result['anchor_strategy'],'median_mapped_seed_midpoint')
        self.assertEqual(len(result['intervals']),7)
        expected=[4000+i*148 for i in range(6)]
        for (a,b),p in zip(result['intervals'][1:],expected):
            self.assertLessEqual(a,p+2)
            self.assertGreaterEqual(b,p+26)

    def test_repetitive_DNA_needs_independent_family_and_context_support(self):
        base=dict(rt_family=True,partner_family=True,compatible_context=True,
            array_status='accepted_array_candidate',copies=6)
        self.assertEqual(classify_locus(**base)[0],'supported_ART_locus_candidate')
        for evidence in ('rt_family','partner_family','compatible_context'):
            altered=dict(base);altered[evidence]=False
            self.assertEqual(classify_locus(**altered)[0],'ART_locus_not_supported')

    def test_all_gap_chance_seed_is_excluded_without_losing_missing_copy_array(self):
        record=SeqIO.read(ROOT/'results/art_validation/loci/missing_copy_2.gb','genbank')
        cfg=json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp') as temp:
            out=Path(temp);(out/'alignments').mkdir()
            result=align_boundaries(str(record.seq[:6000]),[3999,4207,4323,4623,4831,5039],10,cfg,out,'all_gap_outlier')
        self.assertEqual(result['status'],'aligned')
        self.assertEqual(result['excluded_empty_alignment_copies'],[3])
        self.assertEqual(len(result['intervals']),5)

    def test_incomplete_upstream_is_not_complete_locus_even_with_supported_array(self):
        for copies,ast in ((0,'insufficient_recovered_copies'),(6,'accepted_array_candidate')):
            status,reasons=classify_locus(rt_family=True,partner_family=True,compatible_context=True,
                array_status=ast,copies=copies,upstream_shortfall=5200)
            self.assertEqual(status,'incomplete_ART_locus_assessment')
            self.assertIn('upstream_sequence_incomplete',reasons)

    def test_two_copy_or_review_array_never_silently_becomes_complete_locus(self):
        for copies,ast in ((2,'accepted_array_candidate'),(5,'statistically_supported_requires_review')):
            self.assertNotEqual(classify_locus(rt_family=True,partner_family=True,compatible_context=True,
                array_status=ast,copies=copies)[0],'supported_ART_locus_candidate')


if __name__=='__main__':unittest.main()
