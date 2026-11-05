import io
import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import zstandard
from art_type1_public_search import locate_peptide, accession_from_contig_name, normalize_hmmer
from art_type1_archive_leads import scan_prefix
from art_type1_pipeline import best_hits


class TestPublicMappings(unittest.TestCase):
    def test_all_frames_and_reverse(self):
        coding=Seq("ATG"+"GCTGACTGGAAG"*5)
        peptide=str(coding.translate())
        for frame in range(3):
            seq=Seq("C"*frame)+coding+Seq("TAAG")
            for strand in (1,-1):
                source=seq if strand==1 else seq.reverse_complement()
                maps=locate_peptide(SeqRecord(source),peptide)
                self.assertEqual(len(maps),1)
                m=maps[0];got=source[m["start0"]:m["end0"]]
                if m["strand"]==-1:got=got.reverse_complement()
                self.assertEqual(got,coding)
                self.assertEqual(m["strand"],strand)

    def test_duplicate_mapping_is_not_arbitrarily_chosen(self):
        dna=Seq("ATGGCTGACTGGAAGTAG"+"ATGGCTGACTGGAAGTAG")
        self.assertEqual(len(locate_peptide(SeqRecord(dna),"MADWK")),2)

    def test_alternative_initiator_and_code15(self):
        rec=SeqRecord(Seq("GTGGCTGACTGGAAGTAA"))
        self.assertTrue(locate_peptide(rec,"MADWK")[0]["initiator_normalized"])
        rec=SeqRecord(Seq("ATGGCTTAGTGGAAGTAA"))
        m=locate_peptide(rec,"MAQWK")
        self.assertEqual(len(m),1);self.assertEqual(m[0]["genetic_code"],15)
        self.assertFalse(locate_peptide(rec,"MAVWK"))

    def test_real_accession_required(self):
        self.assertEqual(accession_from_contig_name("x-OSOD01046483.1-y"),"OSOD01046483.1")
        self.assertIsNone(accession_from_contig_name("ERZ2641806.86530-NODE-86530"))

    def test_complete_fasta_in_partial_compressed_stream(self):
        content=b">run_1\nACGT\n>run_2\nTGGCAA\n"
        raw=zstandard.ZstdCompressor().compress(content)
        self.assertEqual(str(scan_prefix(raw,"run_1").seq),"ACGT")
        self.assertEqual(str(scan_prefix(raw,"run_2").seq),"TGGCAA")
        self.assertIsNone(scan_prefix(raw[:-5],"run_2"))
        self.assertIsNone(scan_prefix(raw,"missing"))

    def test_rvt_first_domain_defines_ntd(self):
        h=dict(target="x",family="RVT_1",domain_evalue=1e-7,domain_bits=80,protein_envelope_start0=200)
        earlier=dict(h,domain_bits=40,protein_envelope_start0=100)
        self.assertEqual(best_hits([h,earlier],{})["RVT_1"]["x"]["protein_envelope_start0"],100)

    def test_hmmer_missing_identifier_fails(self):
        with self.assertRaises(ValueError):
            normalize_hmmer({"result":{"hits":[{"name":"unknown"}]}},{"query_id":"q","family":"ART_RT"})

    def test_hmmer_real_metadata_accession(self):
        rows=normalize_hmmer({"result":{"hits":[{"name":"031054442","acc":None,
            "metadata":{"accession":"MGYP003692711065"},"evalue":1e-50,"score":210}]}},
            {"query_id":"q","family":"Type_I_partner"})
        self.assertEqual(rows[0]["protein_accession"],"MGYP003692711065")


if __name__=="__main__":unittest.main()
