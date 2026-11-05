"""Coordinate and annotation tests with synthetic loci, not HMM rediscovery."""
import io
import argparse
from contextlib import redirect_stdout
import csv
import hashlib
import json
import random
from pathlib import Path
import tempfile
import unittest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation,
                            SeqFeature, SimpleLocation)
from Bio.SeqRecord import SeqRecord

from extract_art_loci import (Segment, complement_intervals, extract_record,
                             plan_locus, run, source_pieces, write_gff)


def record(length=20000, topology='linear'):
    rng = random.Random(123)
    result = SeqRecord(Seq(''.join(rng.choices('ACGT', k=length))), id='TEST.1')
    result.annotations = dict(molecule_type='DNA', topology=topology)
    return result


def gene(start, end, strand=1, name='RT'):
    return SeqFeature(SimpleLocation(start, end, strand=strand), type='CDS',
        qualifiers=dict(protein_id=[name], codon_start=['1'], translation=['ORIGINAL']))


class Coordinates(unittest.TestCase):
    def test_random_interval_roundtrips_both_strands(self):
        rng = random.Random(71)
        for strand in (-1, 1):
            seg = Segment(213, 1500, 83, 1370, strand)
            for i in range(200):
                a = rng.randrange(213, 1500)
                b = rng.randrange(a+1, 1501)
                self.assertEqual(seg.locus_interval(*seg.source_interval(a, b)), (a, b))
            with self.assertRaises(ValueError):
                seg.source_interval(212, 215)

    def test_reverse_upstream_is_higher_original_coordinates(self):
        rec = record()
        rt = gene(8000, 9500, -1)
        rec.features = [rt, gene(7000, 8020, -1, 'overlapping_partner'), gene(10000, 10450, 1, 'opposite')]
        plan = plan_locus(rec, rt, [rec.features[1]])
        loc = extract_record(rec, rt, plan, 'REVERSE')
        self.assertEqual(str(loc.seq[:6000]), str(rec.seq[9500:15500].reverse_complement()))
        self.assertEqual(plan['segments'], [Segment(2000, 15500, 0, 13500, -1)])
        mapped = {f.qualifiers['protein_id'][0]: f for f in loc.features}
        self.assertEqual(mapped['RT'].location, SimpleLocation(6000, 7500, strand=1))
        self.assertEqual(mapped['overlapping_partner'].extract(loc.seq), rec.features[1].extract(rec.seq))
        self.assertEqual(mapped['opposite'].location.strand, -1)
        self.assertEqual(rec.features[0].qualifiers['translation'], ['ORIGINAL'])

    def test_linear_boundary_clips_without_wrapping(self):
        rec = record(12000)
        rt = gene(100, 1600, 1)
        rec.features = [rt]
        plan = plan_locus(rec, rt)
        self.assertFalse(plan['wrapped'])
        self.assertEqual(plan['upstream_available_nt'], 100)
        self.assertEqual(plan['upstream_shortfall_nt'], 5900)
        self.assertTrue(plan['source_origin_reached'])
        self.assertEqual(plan['segments'], [Segment(0, 7600, 0, 7600, 1)])

    def test_reverse_at_right_boundary_clips_upstream(self):
        rec = record(12000)
        rt = gene(10300, 11800, -1)
        rec.features = [rt]
        plan = plan_locus(rec, rt)
        loc = extract_record(rec, rt, plan, 'EDGE')
        self.assertEqual(plan['upstream_available_nt'], 200)
        self.assertEqual(plan['upstream_shortfall_nt'], 5800)
        self.assertEqual(loc.seq[:200], rec.seq[11800:].reverse_complement())
        self.assertEqual(loc.features[0].location.strand, 1)

    def test_downstream_extended_to_whole_partner(self):
        rec = record(30000)
        rt = gene(10000, 11500)
        partner = gene(17200, 20500, name='long_partner')
        rec.features = [rt, partner]
        plan = plan_locus(rec, rt, [partner])
        loc = extract_record(rec, rt, plan, 'EXTENDED')
        self.assertEqual(plan['downstream_available_nt'], 9000)
        self.assertEqual(loc.features[1].extract(loc.seq), partner.extract(rec.seq))
        self.assertEqual(loc.features[1].qualifiers['locus_clipped'], ['false'])

    def test_circular_wrap_requires_explicit_source_topology(self):
        for strand in (-1, 1):
            rec = record(20000, 'circular')
            rt = gene(200, 1700, strand)
            rec.features = [rt]
            plan = plan_locus(rec, rt)
            loc = extract_record(rec, rt, plan, 'CIRCULAR')
            self.assertTrue(plan['wrapped'])
            self.assertEqual(len(loc), 13500)
            self.assertEqual(loc.features[0].extract(loc.seq), rt.extract(rec.seq))
            self.assertEqual(sum(s.source_end0-s.source_start0 for s in plan['segments']), len(loc))
            for seg in plan['segments']:
                seq = rec.seq[seg.source_start0:seg.source_end0]
                if strand == -1:
                    seq = seq.reverse_complement()
                self.assertEqual(seq, loc.seq[seg.locus_start0:seg.locus_end0])
        rec.annotations.pop('topology')
        self.assertFalse(plan_locus(rec, rt)['wrapped'])

    def test_origin_spanning_join_and_unknown_strand_features(self):
        for strand in (-1, 1):
            rec = record(20000, 'circular')
            parts = [SimpleLocation(19800, 20000, strand=strand), SimpleLocation(0, 400, strand=strand)]
            if strand == -1:
                parts.reverse()
            rt = SeqFeature(CompoundLocation(parts), type='CDS', qualifiers={'protein_id': ['joined']})
            rec.features = [rt, SeqFeature(SimpleLocation(500, 900), type='misc_feature')]
            plan = plan_locus(rec, rt, upstream=1200, downstream=1200)
            loc = extract_record(rec, rt, plan, 'JOINED')
            mapped = next(f for f in loc.features if f.type == 'CDS')
            self.assertEqual(mapped.extract(loc.seq), rt.extract(rec.seq))
            self.assertEqual(mapped.location.strand, 1)

    def test_full_contig_cap_never_duplicates_source_bases(self):
        rec = record(3000, 'circular')
        rt = gene(1000, 1600)
        plan = plan_locus(rec, rt)
        self.assertTrue(plan['circular_request_capped_to_one_contig'])
        self.assertEqual(plan['length_nt'], len(rec))
        coordinates = [p for s in plan['segments'] for p in range(s.source_start0, s.source_end0)]
        self.assertEqual(len(coordinates), len(set(coordinates)))
        self.assertEqual(set(coordinates), set(range(len(rec))))

    def test_clipped_CDS_translation_removed_original_retained_and_phase_adjusted(self):
        rec = record()
        rt = gene(8000, 9500)
        clipped = gene(1900, 2500, name='clipped')
        rec.features = [rt, clipped]
        loc = extract_record(rec, rt, plan_locus(rec, rt), 'CLIPPED')
        f = next(f for f in loc.features if f.qualifiers['protein_id'] == ['clipped'])
        self.assertEqual(f.qualifiers['locus_clipped'], ['true'])
        self.assertNotIn('translation', f.qualifiers)
        self.assertEqual(f.qualifiers['codon_start'], ['3'])
        self.assertEqual(clipped.qualifiers['translation'], ['ORIGINAL'])
        self.assertIsInstance(f.location.start, BeforePosition)

    def test_fuzzy_and_joined_annotations_survive_reverse_and_genbank(self):
        rec = record()
        rt = gene(8000, 9500, -1)
        joined = SeqFeature(CompoundLocation([
            SimpleLocation(BeforePosition(10000), 10100, strand=-1),
            SimpleLocation(9800, AfterPosition(9900), strand=-1)]), type='CDS', qualifiers={'protein_id': ['fuzzy']})
        rec.features = [rt, joined]
        loc = extract_record(rec, rt, plan_locus(rec, rt), 'FUZZY')
        mapped = next(f for f in loc.features if f.qualifiers['protein_id'] == ['fuzzy'])
        self.assertEqual(mapped.extract(loc.seq), joined.extract(rec.seq))
        buf = io.StringIO()
        SeqIO.write(loc, buf, 'genbank')
        buf.seek(0)
        reread = SeqIO.read(buf, 'genbank')
        self.assertEqual([str(f.location) for f in reread.features], [str(f.location) for f in loc.features])

    def test_gff_explicit_one_based_coordinates_and_CDS_phase(self):
        rec = record()
        rt = gene(8000, 9500, -1)
        rec.features = [rt]
        loc = extract_record(rec, rt, plan_locus(rec, rt), 'GFF')
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'test.gff3'
            write_gff(path, loc)
            row = next(line.split('\t') for line in path.read_text().splitlines() if not line.startswith('#'))
        self.assertEqual(row[3:5], ['6001', '7500'])
        self.assertEqual(row[6:8], ['+', '0'])

    def test_issue_interval_mapping_across_circular_origin(self):
        segments = [Segment(90, 100, 0, 10, -1), Segment(80, 90, 10, 20, -1)]
        pieces = source_pieces(8, 13, segments)
        self.assertEqual([(p['source_start0'], p['source_end0']) for p in pieces], [(90, 92), (87, 90)])
        self.assertEqual(complement_intervals([(3, 8), (6, 10), (12, 14)], 17), [(0, 3), (10, 12), (14, 17)])

    def test_unknown_RT_orientation_is_not_invented(self):
        rec = record()
        rt = gene(8000, 9500, None)
        with self.assertRaisesRegex(ValueError, 'orientation is unknown'):
            plan_locus(rec, rt)

    def test_end_to_end_ambiguous_DNA_explicit_gap_and_reverse_coordinate_map(self):
        rec = record()
        rt = gene(8000, 9503, -1, 'synthetic_RT')
        coding = Seq('ATG'+'AAA'*499+'TAA')
        rt.qualifiers.update(translation=['M'+'K'*499], transl_table=['11'])
        text = str(rec.seq)
        text = text[:8000]+str(coding.reverse_complement())+text[9503:]
        rec.seq = Seq(text[:10000]+'NNR'+text[10003:])
        rec.features = [rt, SeqFeature(SimpleLocation(10000, 10002), type='assembly_gap',
                                      qualifiers={'estimated_length': ['2']})]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root/'source.gb'
            SeqIO.write(rec, source, 'genbank')
            sources, candidates = root/'sources.tsv', root/'candidates.tsv'
            source_row = dict(accession=rec.id, path=str(source), format='genbank',
                file_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                sequence_sha256=hashlib.sha256(str(rec.seq).encode()).hexdigest())
            candidate_row = dict(locus_id='synthetic', accession=rec.id, RT_protein_id='synthetic_RT',
                RT_start_1based=8001, RT_end_1based=9503, RT_strand=-1,
                discovery_routes='ART_protein_first', RT_profile_coverage=1.0)
            for path, row in ((sources, source_row), (candidates, candidate_row)):
                with path.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, list(row), delimiter='\t')
                    writer.writeheader()
                    writer.writerow(row)
            output = root/'out'
            with redirect_stdout(io.StringIO()):
                run(argparse.Namespace(output=output, candidates=candidates, sources=sources,
                                       upstream=6000, downstream=6000))
            summary = json.loads((output/'summary.json').read_text())
            self.assertEqual(summary['ambiguous_nt'], 3)
            self.assertEqual(summary['issues_by_kind']['explicit_source_assembly_gap'], 1)
            with (output/'sequence_and_annotation_issues.tsv').open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            ambiguity = next(r for r in rows if r['kind'] == 'ambiguous_sequence')
            details = json.loads(ambiguity['details_json'])
            self.assertEqual(details['symbols'], ['N', 'Y'])
            self.assertEqual(details['source_pieces'][0]['source_start0'], 10000)
            self.assertEqual(details['source_pieces'][0]['source_end0'], 10003)
            self.assertTrue(json.loads((output/'verification.json').read_text())['passed'])
            self.assertEqual(source.read_bytes(), (output/'original_sources'/'TEST.1.gb').read_bytes())


if __name__ == '__main__':
    unittest.main()
