"""Tests for independent ART repeat searches, boundaries, nulls and coordinates."""
import itertools
import json
from pathlib import Path
import random
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
import numpy as np

from art_array_detector import (ROOT, approximate_positions, best_chain, boundary_estimate,
    conserved_block, context_flags, decode, dense_run, distinct_arrays, dna_array, empirical,
    initial_scan, initial_support, nonoverlap, null_batch, pwm_scores, refinement_scan,
    shuffled, stable_seed, words, align_boundaries)
from art_array_initial_candidates import export_initial_candidates
from extract_art_loci import write_json, write_tsv, read_tsv

CFG = json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']


def planted(positions, length=6000, repeat='GACTCATGTGAATCGCATGTAC', seed=17):
    rng = random.Random(seed)
    data = list(''.join(rng.choices('ACGT', k=length)))
    for i, p in enumerate(positions):
        data[p:p+len(repeat)] = repeat
    return ''.join(data)


class DetectorTests(unittest.TestCase):
    def test_encoded_words_and_mismatches_include_imperfect_copies(self):
        sequence = 'ACGTACGTACNNACGTACGTTCNNACGTACGGTC'
        index = words(sequence, 10)
        code = next(c for c in index if decode(c,10) == 'ACGTACGTAC')
        self.assertEqual(approximate_positions(index,code,10,0), [0])
        self.assertEqual(approximate_positions(index,code,10,1), [0,12])
        self.assertEqual(approximate_positions(index,code,10,2), [0,12,24])
        self.assertFalse(any('N' in sequence[p:p+10] for ps in index.values() for p in ps))

    def test_nonoverlapping_matches(self):
        self.assertEqual(nonoverlap([0,1,9,10,20,20],10), [0,10,20])

    def test_single_skipped_copy_and_regular_spacing(self):
        c = best_chain([100,300,700,900], 60,600, missing=1)
        self.assertEqual(c['positions'], [100,300,700,900])
        self.assertEqual(c['missing_copies'],1)
        self.assertAlmostEqual(c['base_spacing'],200)
        c = best_chain([0,100,300,500,900], 90,110, missing=1)
        self.assertLess(len(c['positions']),5)

    def test_dynamic_programming_matches_exhaustive_small_cases(self):
        rng = random.Random(6)
        for trial in range(12):
            positions = sorted(rng.sample(range(0,1300),6))
            best = 1
            for count in range(2,7):
                for chain in itertools.combinations(positions,count):
                    gaps = np.diff(chain)
                    for skip in [-1]+list(range(len(gaps))):
                        factors = np.ones(len(gaps))
                        if skip >= 0:
                            factors[skip] = 2
                        lo = max(60, max(gaps/(factors*1.3)))
                        hi = min(600, min(gaps/(factors*.7)))
                        if lo <= hi+1e-8:
                            best = max(best,count)
            self.assertEqual(len(best_chain(positions,60,600,missing=1)['positions']),best)

    def test_dense_repeat_veto_and_initial_seed_filters(self):
        motif = 'GACTCATGTGAATCGCATGTAC'
        s = planted(list(range(200,850,50)), length=3000, repeat=motif)
        result = initial_scan(s,CFG)
        self.assertTrue(result['dense_veto'])
        self.assertGreater(result['maximum_dense_run_any_seed'],result['R'])
        low_complexity = initial_scan('A'*3000,CFG)
        self.assertEqual(low_complexity['seeds_tested'],0)

    def test_extended_executes_when_initial_has_no_array(self):
        s = planted([200,402,604,806,1008])
        first = initial_scan(s[-3000:],CFG)
        extended = refinement_scan(s,CFG)
        self.assertLess(first['R'],3)
        self.assertTrue(any(len(c['positions']) >= 5 and c['positions'][-1] < 1500
                            for c in extended['candidates']))
        self.assertFalse(initial_support(first,False,3000,[[200,222],[402,424],[604,626]],14))

    def test_initial_success_at_other_array_does_not_change_extended_only_label(self):
        initial = {'positions':[100,300,500]}
        self.assertFalse(initial_support(initial,True,3000,[[200,222],[402,424],[604,626]],14))
        self.assertTrue(initial_support(initial,True,3000,[[3100,3122],[3300,3322],[3500,3522]],14))

    def test_boundaries_use_full_repeat_not_seed_for_spacer_lengths(self):
        motif = 'GACTCATGTGAATCGCATGTAC'
        starts = [600,802,1004,1206,1408]
        s = planted(starts, repeat=motif)
        seed_positions = [p+6 for p in starts]
        boundary = boundary_estimate(s,seed_positions,10,CFG)
        self.assertGreaterEqual(boundary['length_nt'],len(motif))
        intervals = boundary['intervals']
        for (a,b),(c,d) in zip(intervals,intervals[1:]):
            self.assertEqual(c-a,202)
            self.assertEqual(c-b,202-(b-a))
            self.assertLess(c-b,202-10)

    def test_observed_lengths_are_not_hard_limits(self):
        motif = 'GACTCATGTGAATCGCATGTACGCTAGCGTACGATTCGATCGACTGACTAGCTGATCGTA'
        s = planted([600,1000,1400,1800,2200],repeat=motif)
        b = boundary_estimate(s,[p+20 for p in [600,1000,1400,1800,2200]],10,CFG)
        self.assertGreater(b['length_nt'],49)
        self.assertFalse(CFG['enforce_observed_ranges_as_filters'])
        self.assertTrue(refinement_scan(s,CFG)['candidates'])

    def test_local_block_shuffle_preserves_each_block_composition(self):
        s = 'A'*30+'C'*20+'G'*15+'T'*35+'ACGTACG'
        result = shuffled(s,np.random.default_rng(9),50)
        from collections import Counter
        for i in range(0,len(s),50):
            self.assertEqual(Counter(result[i:i+50]),Counter(s[i:i+50]))
        self.assertEqual(result,shuffled(s,np.random.default_rng(9),50))

    def test_shuffles_rescan_seeds_and_reproduce_with_same_seed(self):
        cfg = dict(CFG, initial_top_seeds=4)
        s = planted([100,302,504],length=900)
        task = (s,'initial',cfg,[stable_seed('unit',i) for i in range(2)])
        a = null_batch(task)
        self.assertEqual(a,null_batch(task))
        for seed,score in a:
            self.assertEqual(score,initial_scan(shuffled(s,np.random.default_rng(seed)),cfg)['R'])

    def test_strict_null_max_comparison_and_finite_empirical_p(self):
        t = empirical(10,[8,9,10])
        self.assertFalse(t['exceeds_all_shuffles'])
        self.assertEqual(t['empirical_p'],.5)
        self.assertEqual(empirical(11,[8,9,10])['empirical_p'],.25)

    def test_coding_repeat_and_intervening_gene_filters(self):
        coding = SeqFeature(SimpleLocation(200,1700,strand=-1), type='CDS', qualifiers={'protein_id':['coding']})
        intervals = [[300,322],[501,523],[702,724]]
        flags = context_flags(intervals,[coding],2000,CFG)
        self.assertTrue(flags['coding_repeat'])
        flags = context_flags([[300,322],[502,524],[704,726]],[coding],2000,CFG)
        self.assertFalse(flags['coding_repeat'])
        intervening = SeqFeature(SimpleLocation(1000,1400,strand=1),type='CDS')
        self.assertFalse(context_flags(intervals,[intervening],2000,CFG)['adjacent_to_RT'])

    def test_PWM_does_not_treat_ambiguous_base_as_match(self):
        pwm = np.array([[2,-1,-1,-1],[2,-1,-1,-1]])
        scores = pwm_scores('AANAA',pwm)
        self.assertEqual(scores[0],4)
        self.assertTrue(np.isneginf(scores[1]))
        self.assertTrue(np.isneginf(scores[2]))
        self.assertEqual(scores[3],4)

    def test_one_low_conservation_column_is_allowed(self):
        matrix = np.array([[0,1,2,3,0],[0,1,0,3,0],[0,1,1,3,0]])
        b = conserved_block(matrix,2,.8,1)
        self.assertEqual((b['start'],b['end']),(0,5))
        self.assertIsNone(conserved_block(matrix,2,.8,0))

    def test_MAFFT_indel_repeat_boundaries_map_back_to_input_DNA(self):
        starts = [600,802,1004,1206,1408]
        motif = 'GACTCATGTGAATCGCATGTACGCTAGCGTACGATTCGATCGACTGACTAGCTGATCGTA'
        s = planted(starts,repeat=motif)
        data = list(s)
        insertion = motif[:35]+'C'+motif[35:]
        deletion = motif[:38]+motif[39:]+'G'
        data[starts[1]:starts[1]+len(motif)+1] = insertion
        data[starts[3]:starts[3]+len(motif)] = deletion
        s = ''.join(data)
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)
            (output/'alignments').mkdir()
            result = align_boundaries(s,[p+6 for p in starts],10,CFG,output,'indels')
            self.assertEqual(result['status'],'aligned')
            for (a,b),text in zip(result['intervals'],result['aligned_repeat_sequences']):
                self.assertEqual(s[a:b],text.replace('-',''))
            self.assertTrue(all(n >= len(motif)-1 for n in result['lengths_nt']))
            self.assertTrue(result['gaps_in_repeat'])
            self.assertGreater(len(set(result['lengths_nt'])),1)

    def test_initial_export_retains_dropped_copies_and_reverse_source_coordinates(self):
        starts = [4000,4202,4404,4606,4808]
        seq = planted(starts,length=13503)
        rec = SeqRecord(Seq(seq),id='SYNTHETIC')
        rec.annotations['molecule_type'] = 'DNA'
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inputs, output = root/'inputs',root/'output'
            for folder in (inputs,output/'alignments',output/'details',output/'code'):
                folder.mkdir(parents=True,exist_ok=True)
            SeqIO.write(rec,inputs/'locus.gb','genbank')
            write_tsv(inputs/'locus_manifest.tsv',[dict(locus_id=rec.id,accession='SOURCE',
                genbank_path='locus.gb',RT_locus_start0=6000)],
                ['locus_id','accession','genbank_path','RT_locus_start0'])
            write_tsv(inputs/'coordinate_map.tsv',[dict(locus_id=rec.id,source_start0=50000,
                source_end0=63503,locus_start0=0,locus_end0=13503,orientation=-1)],
                ['locus_id','source_start0','source_end0','locus_start0','locus_end0','orientation'])
            write_tsv(output/'locus_calls.tsv',[dict(locus_id=rec.id,initial_pass=True)],['locus_id','initial_pass'])
            write_tsv(output/'arrays.tsv',[dict(array_id='partial',locus_id=rec.id,copy_count=3,
                repeat_lengths_nt_json='[22,22,22]',array_start0=4000,array_end0=4426)],
                ['array_id','locus_id','copy_count','repeat_lengths_nt_json','array_start0','array_end0'])
            rows = [dict(array_id='partial',locus_start0=p,locus_end0=p+22) for p in starts[:3]]
            write_tsv(output/'repeats.tsv',rows,list(rows[0]))
            write_json(output/'summary.json',{})
            write_json(output/'verification.json',{'passed':True})
            (output/'REPORT.md').write_text('Test report\n')
            write_json(output/'details'/f'{rec.id}.scan.json',{'initial':{'positions':[p+4-3000 for p in starts]}})
            export_initial_candidates(output,inputs,CFG)
            initial = read_tsv(output/'initial_arrays.tsv')[0]
            self.assertEqual(initial['copy_count'],'5')
            self.assertIn('initial_vs_refined_copy_count_disagreement',initial['review_flags'])
            self.assertEqual(initial['matching_delimited_copy_counts'],'3')
            for row in read_tsv(output/'initial_repeats.tsv'):
                piece = json.loads(row['source_segments_json'])[0]
                self.assertEqual(piece['orientation'],-1)
                original = rec.seq.reverse_complement()[piece['source_start0']-50000:piece['source_end0']-50000]
                self.assertEqual(str(original.reverse_complement()),row['sequence'])


if __name__ == '__main__':
    unittest.main()
