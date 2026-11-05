"""Meaningful context/route tests independent of the biological inference results."""
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from art_discovery import ROOT, FAMILIES, Pipeline, neighbors, window_bounds, repeat_hint


def gene(pid, start, end, strand=1):
    return dict(gene_key=f'fixture|{pid}',protein_id=pid,accession='fixture',
                start0=start,end0=end,strand=strand,length_aa=(end-start)//3,
                product='test fixture',translation_qc='exact',partial_location=False,
                source_path='synthetic_test_fixture')


class ContextTests(unittest.TestCase):
    def test_reverse_overlapping_upstream_cds_is_not_skipped(self):
        partner=gene('partner',1000,1500,-1)
        rt=gene('rt',1480,2800,-1)
        other=gene('opposite_strand',2850,3000,1)
        rows=neighbors(partner,[partner,rt,other],'upstream')
        self.assertEqual([g['protein_id'] for g in rows],['rt','opposite_strand'])

    def test_other_strand_genes_count_toward_first_three(self):
        partner=gene('partner',9000,9500)
        genes=[partner,gene('a',8500,8800,-1),gene('b',8000,8300),
               gene('c',7500,7800,-1),gene('fourth',7000,7400)]
        self.assertEqual([g['protein_id'] for g in neighbors(partner,genes,'upstream')],['a','b','c'])

    def test_window_contains_whole_gene_and_clips_at_contig_end(self):
        self.assertEqual(window_bounds(gene('p',100,400),10000,6000),(0,6400))
        self.assertEqual(window_bounds(gene('p',9500,9900),10000,6000),(3500,10000))

    def test_unknown_orientation_has_no_invented_order(self):
        self.assertEqual(neighbors(gene('p',900,1200,0),[gene('a',100,800)],'upstream'),[])

    def test_repeat_hint_on_reverse_strand_with_skipped_copy(self):
        # Exact-word copies separated by 170, 340 and 170 nt. This only tests
        # the hint implementation; it is not a simulated biological dataset.
        word='ATGCTACGTA';sequence=list('N'*1500)
        for start in [100,270,610,780]:sequence[start:start+10]=word
        upstream=Seq(''.join(sequence)).reverse_complement()
        record=SeqRecord(Seq('N'*500)+upstream,id='fixture')
        rt=gene('rt',100,500,-1)
        hint=repeat_hint(record,rt,[rt])
        self.assertEqual(hint['repeat_hint_copies'],4)


class IndependentRouteTests(unittest.TestCase):
    def setUp(self):
        (ROOT/'tmp').mkdir(exist_ok=True)
        self.tmp=tempfile.TemporaryDirectory(dir=ROOT/'tmp')
        p=Pipeline.__new__(Pipeline)
        p.out=Path(self.tmp.name);(p.out/'windows').mkdir()
        p.cfg={'flank_nt':6000}
        p.records={'fixture':SeqRecord(Seq('N'*20000),id='fixture')}
        self.rt=gene('rt',4000,5000)
        self.first=gene('first',5001,5500)
        self.second=gene('second',5501,6000)
        p.genes={g['protein_id']:g for g in [self.rt,self.first,self.second]}
        p.by_contig={'fixture':list(p.genes.values())}
        p.proteins={pid:SeqRecord(Seq('M'+'A'*99),id=pid) for pid in p.genes}
        p.accepted=defaultdict(dict);p.raw=defaultdict(dict)
        p.generic={'rt':{'family':'RVT-Retrons'}}
        p.train={f:[] for f in FAMILIES}
        p.registry=[{'family':f,'status':'provisional'} for f in FAMILIES]
        p.accepted['Type_III_first_partner']['first']={'sequence_bitscore':100}
        p.accepted['Type_III_second_partner']['second']={'sequence_bitscore':100}
        self.pipeline=p

    def tearDown(self):
        self.tmp.cleanup()

    def test_partner_first_recovers_locus_without_ART_profile_hit(self):
        rows=self.pipeline.assemble()
        self.assertEqual(len(rows),1)
        row=rows[0]
        self.assertEqual(row['RT_protein_id'],'rt')
        self.assertTrue(row['found_by_partner_first'])
        self.assertFalse(row['found_by_ART_protein_first'])
        self.assertEqual(row['discovery_routes'],'partner_first')
        self.assertIn('generic_RT_neighbor_without_screened_ART_match',row['notes'])

    def test_type_III_genes_remain_separate_families_at_shared_locus(self):
        row=self.pipeline.assemble()[0]
        self.assertEqual(row['Type_III_first_partner_ids'],'first')
        self.assertEqual(row['Type_III_second_partner_ids'],'second')
        self.assertEqual(row['partner_family_profiles'],'Type_III_first_partner;Type_III_second_partner')

    def test_wrong_orientation_does_not_create_RT_partner_association(self):
        self.pipeline.genes['rt']['strand']=-1
        rows=self.pipeline.assemble()
        self.assertTrue(all(not r['RT_protein_id'] for r in rows))
        self.assertTrue(all(r['candidate_status']=='partner_profile_hit_context_unresolved' for r in rows))

    def test_partner_hit_is_retained_when_no_upstream_RT_is_detected(self):
        self.pipeline.generic={}
        rows=self.pipeline.assemble()
        self.assertEqual(len(rows),2)
        self.assertTrue(all(r['found_by_partner_first'] for r in rows))
        self.assertTrue(all(not r['found_by_ART_protein_first'] for r in rows))

    def test_GNAT_annotation_alone_does_not_assign_type_I_partner(self):
        self.pipeline.accepted=defaultdict(dict)
        self.pipeline.accepted['ART_RT']['rt']={'sequence_bitscore':100,'profile_coverage':1.0}
        self.pipeline.genes['first']['product']='GNAT family acetyltransferase'
        row=self.pipeline.assemble()[0]
        self.assertTrue(row['found_by_ART_protein_first'])
        self.assertFalse(row['found_by_partner_first'])
        self.assertEqual(row['Type_I_partner_ids'],'')


if __name__=='__main__':unittest.main()
