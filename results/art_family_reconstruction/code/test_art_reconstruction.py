import csv
import json
from pathlib import Path
import tempfile
import subprocess
import sys
import unittest
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from art_reconstruction_common import ROOT,Runner,connected_component,longest_contig_representatives,stage_gate
from art_staged_search import load_assembly,find_candidates
from reconstruct_art_family import add_curated_references,broad_class
from extract_art_loci import read_tsv


class ReconstructionTests(unittest.TestCase):
    def test_failed_external_tool_is_not_recorded_as_success(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp') as tmp:
            out=Path(tmp);run=Runner(out)
            with self.assertRaises(subprocess.CalledProcessError):
                run.run([sys.executable,'-c','raise SystemExit(7)'],out/'failure.log')
            record=json.loads((out/'commands.jsonl').read_text().strip())
            self.assertEqual(record['returncode'],7)

    def test_sequence_component_follows_indirect_links_without_architecture(self):
        edges=[('MarsHill','distant_ART'),('distant_ART','no_array_candidate'),('other_RT','unrelated')]
        self.assertEqual(connected_component(iter(edges),'MarsHill'),{'MarsHill','distant_ART','no_array_candidate'})
        self.assertEqual(connected_component([], 'MarsHill'),{'MarsHill'})

    def test_representative_is_longest_contig_not_first_or_longest_protein(self):
        members=[('a','a'),('a','b'),('a','c'),('d','d')]
        meta={'a':{'contig_length_nt':300},'b':{'contig_length_nt':900},'c':{'contig_length_nt':900},'d':{'contig_length_nt':1}}
        self.assertEqual(longest_contig_representatives(members,meta),{'a':'b','d':'d'})

    def test_later_stages_need_validation_and_previous_cost_yield_and_budget(self):
        cfg={'stage3_enabled':False}
        with patch('art_reconstruction_common.validation_gate',return_value={'passed':False,'reasons':['bad_repeat_copies']}):
            self.assertTrue(stage_gate(1,cfg)['allowed'])
            self.assertFalse(stage_gate(2,cfg,{'records_searched':3,'search_seconds':1.,'new_RT_candidates':0})['allowed'])
        with patch('art_reconstruction_common.validation_gate',return_value={'passed':True,'reasons':[]}):
            self.assertFalse(stage_gate(2,cfg)['allowed'])
            self.assertTrue(stage_gate(2,cfg,{'records_searched':3,'search_seconds':1.,'new_RT_candidates':0})['allowed'])
            self.assertFalse(stage_gate(3,cfg,{'records_searched':2,'search_seconds':1.,'new_RT_candidates':0})['allowed'])

    def test_partner_first_keeps_generic_RT_without_strong_ART_and_three_orientations(self):
        genes={k:dict(protein_key=k,gene_key=k,accession='test.1',protein_id=k,start0=a,end0=b,strand=s) for k,a,b,s in [
            ('rt',100,1000,1),('opposed',1010,1100,-1),('partner',1120,1900,1)]}
        cfg=dict(sequence_evalue=1e-5,domain_evalue=1e-5,minimum_ART_profile_coverage=.6)
        base=dict(sequence_evalue=1e-20,domain_evalue=1e-20,profile_coverage=.9,domain_bits=50,protein_envelope_start0=10)
        hits=[dict(base,target='rt',search_family='RVT_1'),dict(base,target='partner',search_family='Type_II_partner')]
        routes,art,rvt,p,weak,up,orphan=find_candidates(genes,{'test.1':list(genes.values())},hits,cfg)
        self.assertEqual(routes['rt'],{'partner_first'});self.assertFalse(art);self.assertFalse(orphan)
        self.assertEqual(len(up),2);self.assertFalse(up[0]['same_orientation'])

    def test_protein_only_reference_needs_curation_and_never_gets_locus_id(self):
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp') as tmp:
            folder=Path(tmp);fasta=folder/'reference.faa';manifest=folder/'references.tsv'
            SeqIO.write(SeqRecord(Seq('M'+'ACDEFGHIKLMNPQRSTVWY'*15),id='NART_test',description=''),fasta,'fasta')
            manifest.write_text('protein_id\tbroad_class\tcuration_evidence\tsource_url\nNART_test\tNART\t\thttps://example.org\n')
            with self.assertRaises(ValueError):add_curated_references({}, {},fasta,manifest)
            manifest.write_text('protein_id\tbroad_class\tcuration_evidence\tsource_url\nNART_test\tNART\ttest fixture label only\thttps://example.org\n')
            pool={};meta={};added=add_curated_references(pool,meta,fasta,manifest)
            self.assertEqual(added[0]['broad_class'],'NART');self.assertEqual(added[0]['locus_id'],'')

    def test_viral_calls_preserve_source_and_use_one_coordinate_system(self):
        manifest=read_tsv(ROOT/'config/art_stage2_assemblies.tsv');row=next(r for r in manifest if r['annotation_method']=='pyrodigal_gv')
        original=SeqIO.read(ROOT/row['annotation_path'],'genbank')
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp') as tmp:
            rec,meta=load_assembly(row,Path(tmp));source_count=len(original.features)
            self.assertEqual(len(rec.seq),len(original.seq));self.assertEqual(str(rec.seq),str(original.seq))
            self.assertEqual(rec.features[:source_count],original.features)
            cds=[f for f in rec.features if f.type=='CDS'];self.assertGreater(len(cds),100)
            for f in cds:
                self.assertTrue(0<=int(f.location.start)<int(f.location.end)<=len(rec))
                if type(f.location.start).__name__=='ExactPosition' and type(f.location.end).__name__=='ExactPosition':
                    translated=str(f.extract(rec.seq).translate(table=int(f.qualifiers['transl_table'][0]),cds=True))
                    self.assertEqual(translated,f.qualifiers['translation'][0])

    def test_annotation_nucleotide_mismatch_is_rejected(self):
        row=read_tsv(ROOT/'config/art_stage2_assemblies.tsv')[0]
        with tempfile.TemporaryDirectory(dir=ROOT/'tmp') as tmp:
            folder=Path(tmp);fasta=folder/'wrong.fna';SeqIO.write(SeqRecord(Seq('ACGT'),id=row['contig_id'],description=''),fasta,'fasta')
            with self.assertRaisesRegex(ValueError,'sequences differ'):
                load_assembly(dict(row,nucleotide_path=fasta.relative_to(ROOT).as_posix()),folder)


if __name__=='__main__':unittest.main()
