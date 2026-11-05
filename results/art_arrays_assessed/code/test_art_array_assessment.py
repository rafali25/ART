import itertools
import json
from collections import Counter
from pathlib import Path
import random
import unittest
from unittest.mock import patch

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,SimpleLocation

from art_array_assessment import (ROOT,dinucleotide_shuffle,dinucleotide_graph,sensitivity_batch,
    bonferroni,main_status,quality_flags,short_period_scan,aligned_pwm)
from art_array_detector import refinement_scan,shuffled,pwm_scores

CFG = json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
ASSESS = json.loads((ROOT/'config/art_array_assessment.json').read_text())


def pairs(text):
    return Counter(zip(text,text[1:]))


class AssessmentTests(unittest.TestCase):
    def test_dinucleotide_shuffle_preserves_counts_endpoints_and_ambiguity(self):
        rng = random.Random(18)
        text = ''.join(rng.choices('ACGT',k=1200))+'NNRAAACCTGATCGTTAACGNN'+'AGT'*30
        original_pairs = pairs(text)
        for seed in range(20):
            result = dinucleotide_shuffle(text,np.random.default_rng(seed))
            self.assertEqual(pairs(result),original_pairs)
            self.assertEqual(Counter(result),Counter(text))
            self.assertEqual((result[0],result[-1]),(text[0],text[-1]))
            self.assertEqual([i for i,b in enumerate(result) if b not in 'ACGT'],
                             [i for i,b in enumerate(text) if b not in 'ACGT'])
        self.assertEqual(dinucleotide_shuffle(text,np.random.default_rng(7)),
                         dinucleotide_shuffle(text,np.random.default_rng(7)))

    def test_dinucleotide_uniformity_against_exhaustive_short_strings(self):
        text = 'AAGACAGC'
        possibilities = {''.join(p) for p in itertools.permutations(text) if p[0] == text[0] and p[-1] == text[-1]
                         and pairs(''.join(p)) == pairs(text)}
        counts = Counter(dinucleotide_shuffle(text,np.random.default_rng(i)) for i in range(6000))
        self.assertEqual(set(counts),possibilities)
        expected = 6000/len(possibilities)
        self.assertTrue(all(abs(n-expected) < .12*expected for n in counts.values()))

    def test_constant_or_unique_euler_trail_is_not_forced_to_change(self):
        for text in ('AAAAAA','ACGT','A','NNNN','AANNCC'):
            result = dinucleotide_shuffle(text,np.random.default_rng(1))
            self.assertEqual(pairs(text),pairs(result))

    def test_null_batch_rediscovers_seeds_in_each_background(self):
        rng = random.Random(42)
        seq = ''.join(rng.choices('ACGT',k=900))
        for kind in ('local_blocks','dinucleotide'):
            rows = sensitivity_batch((seq,kind,CFG,[125,126]))
            for seed,score in rows:
                r = np.random.default_rng(seed)
                permuted = dinucleotide_shuffle(seq,r) if kind == 'dinucleotide' else shuffled(seq,r,50)
                self.assertEqual(score,refinement_scan(permuted,CFG)['score'])

    def test_distinct_locus_motifs_are_discovered_without_a_universal_repeat(self):
        motifs = ['GACTCATGTGAATCGCATGTAC','CCGTTAGCACCTGAGTTCGACT']
        seeds = []
        for motif in motifs:
            rng = random.Random(17)
            data = list(''.join(rng.choices('ACGT',k=2000)))
            for p in [200,402,604,806,1008]:
                data[p:p+len(motif)] = motif
            result = refinement_scan(''.join(data),CFG)
            matching = [c for c in result['candidates'] if len(c['positions']) >= 5]
            self.assertTrue(matching)
            seeds.append(matching[0]['seed'])
        self.assertNotEqual(seeds[0],seeds[1])

    def test_bonferroni_includes_all_loci_and_two_copy_never_accepted(self):
        self.assertAlmostEqual(bonferroni(.004,9),.036)
        self.assertEqual(bonferroni(.2,9),1.)
        self.assertEqual(main_status(2,0.,[]),'exploratory_two_copy_not_accepted')
        self.assertEqual(main_status(3,.06,[]),'not_supported_after_locus_correction')
        self.assertEqual(main_status(4,.01,[]),'accepted_array_candidate')
        self.assertEqual(main_status(4,.01,['repeat_overlaps_CDS']),'statistically_supported_requires_review')

    def test_context_flags_coding_gene_low_complexity_ambiguity_and_boundaries(self):
        rec = SeqRecord(Seq('A'*100+'NN'+'ACGT'*1500))
        rec.features = [SeqFeature(SimpleLocation(10,70,strand=-1),type='CDS',qualifiers={'protein_id':['overlap']}),
                        SeqFeature(SimpleLocation(400,800,strand=1),type='CDS',qualifiers={'protein_id':['between']})]
        meta = {'upstream_shortfall_nt':'15','source_origin_reached':'True'}
        flags = quality_flags(rec,[[0,30],[90,120]],1000,meta,CFG,ASSESS,
            {'score':100.,'positions':[2,20,40]},1.)['flags']
        for expected in ('repeat_overlaps_CDS','low_complexity_repeat','ambiguous_bases_in_array',
                         'gene_between_array_and_RT','contig_limited_upstream','locus_reaches_source_boundary',
                         'array_touches_upstream_search_boundary','stronger_short_period_pattern'):
            self.assertIn(expected,flags)

    def test_copy_disagreement_requires_review(self):
        self.assertEqual(main_status(3,.001,['initial_vs_final_copy_count_disagreement']),
                         'statistically_supported_requires_review')

    def test_noncoding_gene_between_array_and_RT_is_flagged(self):
        rec = SeqRecord(Seq('ACGT'*500))
        rec.features = [SeqFeature(SimpleLocation(500,575,strand=-1),type='tRNA',qualifiers={'locus_tag':['tRNA_1']})]
        quality = quality_flags(rec,[[100,122],[302,324]],1000,{},CFG,ASSESS,
                                {'score':0.,'positions':[]},50.)
        self.assertIn('gene_between_array_and_RT',quality['flags'])
        self.assertEqual(quality['intervening_gene_ids'],['tRNA_1'])

    def test_short_period_scanner_selects_a_chain_not_just_word_count(self):
        motif = 'GACTCATGTGAATCGCATGTAC'
        sequence = (motif+'TCGACGTA')*8
        result = short_period_scan(sequence,CFG)
        self.assertGreaterEqual(result['copy_count'],6)
        self.assertLess(result['base_spacing'],60)

    def test_repeat_alignment_occurs_before_PWM_rescanning(self):
        motif = 'GACTCATGTGAATCGCATGTAC'
        starts = [200,402,604,806,1008]
        rng = random.Random(17)
        data = list(''.join(rng.choices('ACGT',k=2000)))
        for p in starts: data[p:p+len(motif)] = motif
        sequence = ''.join(data)
        events = []
        def align(*args):
            events.append('alignment')
            return dict(status='aligned',aligned_repeat_sequences=[motif]*5,
                        intervals=[[p,p+len(motif)] for p in starts])
        def score(*args):
            events.append('PWM_scoring')
            return pwm_scores(*args)
        with patch('art_array_assessment.align_boundaries',side_effect=align), \
             patch('art_array_assessment.pwm_scores',side_effect=score):
            model = aligned_pwm(sequence,{'positions':[p+6 for p in starts]},dict(CFG,pwm_shuffles=3),Path('.'),'unit')
        self.assertEqual(events[0],'alignment')
        self.assertEqual(model['status'],'rescanned')
        self.assertEqual(len(model['final_intervals']),5)


if __name__ == '__main__':
    unittest.main()
