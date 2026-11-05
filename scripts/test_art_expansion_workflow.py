"""Meaningful exact-mapping checks for the faster archive prefilter."""
import unittest
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import zstandard
from art_expansion_leads import match_records


class ExactArchivePrefilter(unittest.TestCase):
    peptide='MADWKLSFANGYVAPCEDRHIMQGK'

    def coding(self):
        table=CodonTable.unambiguous_dna_by_id[11].forward_table
        return ''.join(next(c for c,a in table.items() if a==amino) for amino in self.peptide)

    def scan(self,records,complete=True):
        raw=zstandard.ZstdCompressor().compress(records.encode())
        peptide=SeqRecord(Seq(self.peptide),id='lead')
        return match_records(raw,complete,['lead'],{'lead':peptide})[0]

    def test_full_peptide_maps_in_all_frames_and_source_orientations(self):
        for frame in range(3):
            dna='C'*frame+self.coding()+'TAA'
            for strand in [1,-1]:
                seq=dna if strand==1 else str(Seq(dna).reverse_complement())
                hits=self.scan('>source\n'+seq+'\n')
                self.assertEqual(len(hits['lead']),1)
                self.assertEqual(hits['lead'][0][1][0]['strand'],strand)

    def test_genetic_code_15_is_retained_by_DNA_prefilter(self):
        dna=self.coding();i=self.peptide.index('Q')*3
        dna=dna[:i]+'TAG'+dna[i+3:]+'TAA'
        hits=self.scan('>source\n'+dna+'\n')
        self.assertEqual(hits['lead'][0][1][0]['genetic_code'],15)

    def test_unterminated_final_record_is_not_promoted_from_partial_object(self):
        self.assertFalse(self.scan('>last\n'+self.coding()+'\n',complete=False))
        hits=self.scan('>first\n'+self.coding()+'\n>last\nACGT\n',complete=False)
        self.assertEqual(hits['lead'][0][0].id,'first')

    def test_multiple_genomic_copies_remain_ambiguous(self):
        dna=self.coding()+'TAA'+self.coding()+'TAA'
        hits=self.scan('>source\n'+dna+'\n')
        self.assertEqual(len(hits['lead'][0][1]),2)


if __name__=='__main__':unittest.main()
