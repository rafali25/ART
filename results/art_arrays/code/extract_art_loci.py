"""Extract RT-oriented loci without changing source annotations.

All Python/JSON/TSV intervals are zero-based, half-open. GenBank and GFF3
serialize these explicitly as one-based inclusive coordinates.
"""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import quote

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation,
                            ExactPosition, SeqFeature, SimpleLocation)
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
PARTNERS = ('Type_I_partner_ids', 'Type_II_partner_ids',
            'Type_III_first_partner_ids', 'Type_III_second_partner_ids')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=str) + '\n', encoding='utf-8')


def pid(feature):
    return feature.qualifiers.get('protein_id', feature.qualifiers.get('locus_tag', ['']))[0]


@dataclass(frozen=True)
class Segment:
    source_start0: int
    source_end0: int
    locus_start0: int
    locus_end0: int
    orientation: int

    def source_interval(self, start, end):
        if not self.source_start0 <= start <= end <= self.source_end0:
            raise ValueError('Interval is outside mapping segment')
        if self.orientation == 1:
            return (self.locus_start0 + start - self.source_start0,
                    self.locus_start0 + end - self.source_start0)
        return (self.locus_start0 + self.source_end0 - end,
                self.locus_start0 + self.source_end0 - start)

    def locus_interval(self, start, end):
        if not self.locus_start0 <= start <= end <= self.locus_end0:
            raise ValueError('Interval is outside mapping segment')
        if self.orientation == 1:
            return (self.source_start0 + start - self.locus_start0,
                    self.source_start0 + end - self.locus_start0)
        return (self.source_end0 - (end - self.locus_start0),
                self.source_end0 - (start - self.locus_start0))


def enclosing_arc(location, length, circular):
    """Smallest containing arc; linear joins retain their full genomic span."""
    parts = sorted((int(p.start), int(p.end)) for p in location.parts)
    if any(p.ref or p.ref_db for p in location.parts):
        raise ValueError('Remote-reference RT/partner locations require explicit resolution')
    if not parts or any(a < 0 or b > length or a >= b for a, b in parts):
        raise ValueError('Invalid local feature interval')
    if not circular or len(parts) == 1:
        return parts[0][0], max(b for a, b in parts)
    merged = []
    for a, b in parts:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(b, merged[-1][1])
        else:
            merged.append([a, b])
    gaps = [(merged[(i+1) % len(merged)][0] + (length if i == len(merged)-1 else 0) - b, i)
            for i, (a, b) in enumerate(merged)]
    gap, i = max(gaps)
    start = merged[(i+1) % len(merged)][0]
    return start, start + length - gap


def plan_locus(record, rt, partners=(), upstream=6000, downstream=6000):
    if upstream < 0 or downstream < 0:
        raise ValueError('Flank lengths must be nonnegative')
    strand = rt.location.strand
    if strand not in (-1, 1):
        raise ValueError('RT orientation is unknown; cannot orient this locus')
    length = len(record)
    circular = record.annotations.get('topology', '').lower() == 'circular'
    rt_lo, rt_hi = enclosing_arc(rt.location, length, circular)
    core_lo, core_hi = rt_lo, rt_hi
    for partner in partners:
        a, b = enclosing_arc(partner.location, length, circular)
        if circular:
            shift = min((-length, 0, length), key=lambda d: abs((a+b)/2+d-(rt_lo+rt_hi)/2))
            a, b = a+shift, b+shift
        core_lo, core_hi = min(core_lo, a), max(core_hi, b)
    lo = min(core_lo, rt_lo - (upstream if strand == 1 else downstream))
    hi = max(core_hi, rt_hi + (downstream if strand == 1 else upstream))
    wanted_lo, wanted_hi = lo, hi
    if circular:
        if core_hi-core_lo > length:
            raise ValueError('RT and partner core spans more than one contig')
        if hi-lo > length:
            # Do not duplicate genomic bases; retain the complete core first.
            lo = max(core_hi-length, min(lo, core_lo))
            hi = lo+length
    else:
        lo, hi = max(0, lo), min(length, hi)
    plus_segments = []
    cursor = lo
    offset = 0
    while cursor < hi:
        source_start = cursor % length
        count = min(hi-cursor, length-source_start)
        plus_segments.append(Segment(source_start, source_start+count, offset, offset+count, 1))
        offset += count
        cursor += count
    segments = plus_segments if strand == 1 else [
        Segment(s.source_start0, s.source_end0, offset-s.locus_end0, offset-s.locus_start0, -1)
        for s in reversed(plus_segments)]
    upstream_available = rt_lo-lo if strand == 1 else hi-rt_hi
    downstream_available = hi-rt_hi if strand == 1 else rt_lo-lo
    return dict(segments=segments, plus_segments=plus_segments, strand=strand,
                length_nt=offset, upstream_requested_nt=upstream, downstream_requested_nt=downstream,
                upstream_available_nt=upstream_available, downstream_available_nt=downstream_available,
                upstream_shortfall_nt=max(0, upstream-upstream_available),
                downstream_shortfall_nt=max(0, downstream-downstream_available),
                source_topology=record.annotations.get('topology', 'unknown'),
                circularity_support='source_record_topology=circular' if circular else 'none',
                wrapped=(lo < 0 or hi > length) if circular else False,
                source_origin_reached=any(s.source_start0 == 0 for s in segments),
                source_end_reached=any(s.source_end0 == length for s in segments),
                genomic_left_clipped_nt=max(0, -wanted_lo) if not circular else 0,
                genomic_right_clipped_nt=max(0, wanted_hi-length) if not circular else 0,
                circular_request_capped_to_one_contig=circular and wanted_hi-wanted_lo > length)


def feature_snapshot(feature, index):
    loc = feature.location
    return dict(source_feature_index=index, type=feature.type, id=feature.id,
                original_location=str(loc), operator=getattr(loc, 'operator', None),
                parts=[] if loc is None else [dict(start0=int(p.start), end0=int(p.end), strand=p.strand,
                    start_position=repr(p.start), end_position=repr(p.end), ref=p.ref, ref_db=p.ref_db)
                    for p in loc.parts], qualifiers=copy.deepcopy(feature.qualifiers))


def map_feature(feature, index, segments):
    """Map to forward-window coordinates; preserve join order and fuzzy ends.

    Final RT-direction reversal is delegated to Biopython SeqRecord, which also
    flips strands and fuzzy endpoints. Clipped CDS translations are archived in
    the source snapshot rather than mislabelled as translations of a fragment.
    """
    if feature.location is None:
        return None
    if any(p.ref or p.ref_db for p in feature.location.parts):
        raise ValueError('Remote feature location cannot be remapped implicitly')
    mapped = []
    first_offset = None
    consumed = 0
    for part in feature.location.parts:
        matches = []
        for seg in segments:
            a, b = max(int(part.start), seg.source_start0), min(int(part.end), seg.source_end0)
            if a < b:
                matches.append((a, b, seg))
        matches.sort(key=lambda item: item[0], reverse=part.strand == -1)
        for a, b, seg in matches:
            if first_offset is None:
                first_offset = consumed + (int(part.end)-b if part.strand == -1 else a-int(part.start))
            start = part.start if a == int(part.start) else BeforePosition(a)
            end = part.end if b == int(part.end) else AfterPosition(b)
            mapped.append(SimpleLocation(start, end, strand=part.strand) +
                          (seg.locus_start0-seg.source_start0))
        consumed += len(part)
    if not mapped:
        return None
    location = mapped[0] if len(mapped) == 1 else CompoundLocation(mapped, getattr(feature.location, 'operator', 'join'))
    result = copy.deepcopy(feature)
    result.location = location
    result.qualifiers['source_feature_index'] = [str(index)]
    clipped = len(location) != len(feature.location)
    result.qualifiers['locus_clipped'] = [str(clipped).lower()]
    if clipped:
        result.qualifiers.setdefault('note', []).append('Clipped by locus extraction; full annotation in source snapshot.')
        if feature.type == 'CDS':
            result.qualifiers.pop('translation', None)
            original_phase = int(feature.qualifiers.get('codon_start', ['1'])[0])-1
            result.qualifiers['codon_start'] = [str((original_phase-first_offset) % 3+1)]
    return result


def extract_record(record, rt, plan, locus_id):
    seq = Seq('').join([record.seq[s.source_start0:s.source_end0] for s in plan['plus_segments']])
    locus = SeqRecord(seq, id=locus_id, name=locus_id, description=f'{record.id}; RT={pid(rt)}; RT-oriented locus')
    locus.annotations = {k: copy.deepcopy(v) for k, v in record.annotations.items()
                         if k in ('molecule_type', 'organism', 'taxonomy', 'source', 'data_file_division')}
    locus.annotations.update(molecule_type='DNA', topology='linear')
    locus.annotations['comment'] = 'Derived linear locus. Original record and all original annotations are archived separately.'
    for index, feature in enumerate(record.features):
        mapped = map_feature(feature, index, plan['plus_segments'])
        if mapped is not None:
            locus.features.append(mapped)
    if plan['strand'] == -1:
        locus = locus.reverse_complement(id=True, name=True, description=True, features=True, annotations=True)
    # Biopython stores locations internally as 0-based half-open; the writer
    # converts them to GenBank 1-based inclusive locations, including fuzziness.
    return locus


def complement_intervals(intervals, length):
    cursor = 0
    result = []
    for a, b in sorted(intervals):
        if a > cursor:
            result.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < length:
        result.append((cursor, length))
    return result


def source_pieces(start, end, segments):
    result = []
    for seg in segments:
        a, b = max(start, seg.locus_start0), min(end, seg.locus_end0)
        if a < b:
            x, y = seg.locus_interval(a, b)
            result.append(dict(source_start0=x, source_end0=y, orientation=seg.orientation,
                               locus_start0=a, locus_end0=b))
    return result


def write_gff(path, locus):
    """One line per location part. GFF coordinates: start0+1, end0."""
    esc = lambda x: quote(str(x), safe='.:^*$@!+_?-|')
    lines = ['##gff-version 3', f'##sequence-region {locus.id} 1 {len(locus)}',
             '# Coordinates: 1-based inclusive; mapping tables: 0-based half-open.']
    for feature in locus.features:
        index = feature.qualifiers['source_feature_index'][0]
        feature_id = f'{locus.id}.f{index}'
        offset = 0
        phase0 = int(feature.qualifiers.get('codon_start', ['1'])[0])-1
        for part in feature.location.parts:
            attrs = [('ID', [feature_id])]
            attrs += [(k, v) for k, v in feature.qualifiers.items() if k != 'ID']
            if feature.type == 'source':
                attrs.append(('original_feature_type', ['source']))
            text = ';'.join(esc(k)+'='+','.join(esc(v) for v in values) for k, values in attrs)
            phase = str((phase0-offset) % 3) if feature.type == 'CDS' and part.strand in (-1, 1) else '.'
            lines.append('\t'.join((locus.id, 'source_annotation', 'region' if feature.type == 'source' else feature.type,
                str(int(part.start)+1), str(int(part.end)), '.',
                {1: '+', -1: '-'}.get(part.strand, '.'), phase, text)))
            offset += len(part)
    path.write_text('\n'.join(lines)+'\n', encoding='utf-8')


def rt_review(record, rt, candidate):
    reasons = []
    if any(type(p.start) is not ExactPosition or type(p.end) is not ExactPosition for p in rt.location.parts):
        reasons.append('fuzzy_or_partial_source_location')
    if len(rt.location.parts) > 1:
        reasons.append('compound_source_location_review_join_or_split_annotation')
    if any(int(p.start) == 0 or int(p.end) == len(record) for p in rt.location.parts):
        reasons.append('RT_touches_source_contig_boundary')
    if set(str(rt.extract(record.seq)).upper()) - set('ACGT'):
        reasons.append('ambiguous_RT_DNA')
    if any(k in rt.qualifiers for k in ('pseudo', 'pseudogene', 'exception', 'transl_except')):
        reasons.append('source_exception_or_pseudogene_annotation')
    if any(f.type in ('gap', 'assembly_gap') and f.location is not None and
           any(max(int(a.start), int(b.start)) < min(int(a.end), int(b.end))
               for a in rt.location.parts for b in f.location.parts) for f in record.features):
        reasons.append('RT_overlaps_source_gap_annotation')
    table = int(rt.qualifiers.get('transl_table', ['1'])[0])
    try:
        translated = rt.translate(record.seq, cds=True)
        annotated = rt.qualifiers.get('translation', [''])[0].replace(' ', '')
        translation_ok = str(translated) == annotated if annotated else True
        if not translation_ok:
            reasons.append('source_translation_mismatch')
    except Exception as error:
        translation_ok = False
        reasons.append('source_CDS_translation_check:' + str(error))
    coverage = candidate.get('RT_profile_coverage', '')
    if coverage and float(coverage) < 0.9:
        reasons.append('RT_profile_coverage_below_review_threshold_0.9')
    return dict(translation_table=table, source_translation_verified=translation_ok,
                review_required=bool(reasons), review_reasons=reasons,
                interpretation='review_flags_do_not_establish_biological_truncation')


def run(args):
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    for sub in ('loci', 'original_sources', 'code'):
        (out/sub).mkdir(exist_ok=True)
    candidates = read_tsv(args.candidates)
    sources = {r['accession']: r for r in read_tsv(args.sources)}
    cache = {}
    manifests, maps, genes, issues, reviews, verified = [], [], [], [], [], []
    records_out = []
    skipped = []
    for candidate in candidates:
        if not candidate['RT_protein_id']:
            skipped.append(dict(locus_id=candidate['locus_id'], reason='partner_only_locus_without_candidate_RT'))
            continue
        accession, rt_id = candidate['accession'], candidate['RT_protein_id']
        if accession not in cache:
            source = sources[accession]
            path = ROOT/source['path']
            if sha(path.read_bytes()) != source['file_sha256']:
                raise ValueError(f'Source changed: {path}')
            with path.open(encoding='utf-8') as handle:
                record = next(r for r in SeqIO.parse(handle, source['format']) if r.id == accession)
            if sha(str(record.seq).encode()) != source['sequence_sha256']:
                raise ValueError(f'Sequence checksum mismatch: {accession}')
            snapshot_name = accession + path.suffix
            shutil.copyfile(path, out/'original_sources'/snapshot_name)
            write_json(out/'original_sources'/f'{accession}.annotations.json', dict(
                id=record.id, description=record.description, annotations=record.annotations,
                dbxrefs=record.dbxrefs, features=[feature_snapshot(f, i) for i, f in enumerate(record.features)]))
            cache[accession] = record
        record = cache[accession]
        cds = {}
        for feature in record.features:
            if feature.type == 'CDS' and pid(feature):
                if pid(feature) in cds:
                    raise ValueError(f'Duplicate CDS identifier: {accession} {pid(feature)}')
                cds[pid(feature)] = feature
        if rt_id not in cds:
            raise ValueError(f'RT {rt_id} has no source CDS; supply gene-caller annotations before extraction')
        rt = cds[rt_id]
        if (int(candidate['RT_start_1based'])-1, int(candidate['RT_end_1based']), int(candidate['RT_strand'])) != (
                int(rt.location.start), int(rt.location.end), rt.location.strand):
            raise ValueError(f'Candidate coordinates disagree with source annotation: {rt_id}')
        partner_ids = {v for field in PARTNERS for v in candidate.get(field, '').split(';') if v}
        if partner_ids - cds.keys():
            raise ValueError(f'Missing partner CDSs: {partner_ids-cds.keys()}')
        plan = plan_locus(record, rt, [cds[p] for p in sorted(partner_ids)], args.upstream, args.downstream)
        name = f'{accession}__{rt_id}'
        locus = extract_record(record, rt, plan, name)
        records_out.append(locus)
        # Sentinel suffix prevents accession version '.1' being replaced.
        base = out/'loci'/(name + '.locus')
        SeqIO.write(locus, base.with_suffix('.fna'), 'fasta')
        SeqIO.write(locus, base.with_suffix('.gb'), 'genbank')
        write_gff(base.with_suffix('.gff3'), locus)
        # Check every base against its original source segment and round trip.
        for number, seg in enumerate(plan['segments']):
            expected = record.seq[seg.source_start0:seg.source_end0]
            if seg.orientation == -1:
                expected = expected.reverse_complement()
            assert locus.seq[seg.locus_start0:seg.locus_end0] == expected
            assert seg.locus_interval(*seg.source_interval(seg.source_start0, seg.source_end0)) == (
                seg.source_start0, seg.source_end0)
            maps.append(dict(locus_id=name, accession=accession, segment_index=number, **asdict(seg)))
        rt_mapped = next(f for f in locus.features if f.type == 'CDS' and pid(f) == rt_id)
        assert rt_mapped.location.strand == 1
        assert rt_mapped.extract(locus.seq) == rt.extract(record.seq)
        assert rt_mapped.qualifiers['locus_clipped'] == ['false']
        checked_features = 0
        for f in locus.features:
            source_index = int(f.qualifiers['source_feature_index'][0])
            original = record.features[source_index]
            if f.qualifiers['locus_clipped'] == ['false']:
                assert f.extract(locus.seq) == original.extract(record.seq), (name, source_index)
                checked_features += 1
            if f.type in ('CDS', 'gene'):
                genes.append(dict(locus_id=name, accession=accession, source_feature_index=source_index,
                    feature_type=f.type, protein_id=pid(f), role='ART_RT' if pid(f) == rt_id else
                    ('candidate_partner' if pid(f) in partner_ids else 'context_gene'),
                    original_location=str(original.location), source_strand=original.location.strand,
                    locus_start0=int(f.location.start), locus_end0=int(f.location.end),
                    locus_strand=f.location.strand, locus_location=str(f.location),
                    locus_clipped=f.qualifiers['locus_clipped'][0],
                    qualifiers_json=json.dumps(f.qualifiers, sort_keys=True)))
        assert partner_ids <= {pid(f) for f in locus.features if f.qualifiers['locus_clipped'] == ['false']}
        for match in re.finditer('[^ACGT]+', str(locus.seq).upper()):
            a, b = match.span()
            issues.append(dict(locus_id=name, kind='ambiguous_sequence', locus_start0=a, locus_end0=b,
                details_json=json.dumps(dict(symbols=sorted(set(match.group())),
                    source_pieces=source_pieces(a, b, plan['segments'])))))
        for label, types in (('no_CDS_annotation', {'CDS'}),
                             ('no_non_source_annotation', {f.type for f in locus.features if f.type != 'source'})):
            covered = [(int(p.start), int(p.end)) for f in locus.features if f.type in types for p in f.location.parts]
            for a, b in complement_intervals(covered, len(locus)):
                issues.append(dict(locus_id=name, kind=label, locus_start0=a, locus_end0=b,
                    details_json=json.dumps(dict(interpretation='annotation_coverage_only_not_an_assembly_gap',
                        source_pieces=source_pieces(a, b, plan['segments'])))))
        for feature in locus.features:
            if feature.type in ('gap', 'assembly_gap'):
                issues.append(dict(locus_id=name, kind='explicit_source_'+feature.type,
                    locus_start0=int(feature.location.start), locus_end0=int(feature.location.end),
                    details_json=json.dumps(feature.qualifiers)))
            if feature.qualifiers['locus_clipped'] == ['true'] and feature.type != 'source':
                issues.append(dict(locus_id=name, kind='annotation_clipped_by_locus',
                    locus_start0=int(feature.location.start), locus_end0=int(feature.location.end),
                    details_json=json.dumps(dict(source_feature_index=feature.qualifiers['source_feature_index'][0],
                                                feature_type=feature.type, protein_id=pid(feature)))))
        review = rt_review(record, rt, candidate)
        reviews.append(dict(locus_id=name, RT_protein_id=rt_id, **review))
        meta = {k: v for k, v in plan.items() if k not in ('segments', 'plus_segments')}
        meta.update(locus_id=name, accession=accession, RT_protein_id=rt_id,
                    RT_locus_start0=int(rt_mapped.location.start), RT_locus_end0=int(rt_mapped.location.end),
                    partner_ids=';'.join(sorted(partner_ids)), discovery_routes=candidate['discovery_routes'],
                    source_path=sources[accession]['path'], source_file_sha256=sources[accession]['file_sha256'],
                    annotation_source='original_source_record', gene_caller_used=False,
                    ambiguous_nt=sum(len(m.group()) for m in re.finditer('[^ACGT]+', str(locus.seq).upper())),
                    fasta_path=f'loci/{name}.fna', genbank_path=f'loci/{name}.gb', gff3_path=f'loci/{name}.gff3')
        manifests.append(meta)
        write_json(base.with_suffix('.map.json'), dict(coordinate_convention='zero_based_half_open', **meta,
                   segments=[asdict(s) for s in plan['segments']], RT_review=review))
        # Reparse GenBank to verify conversion, joined features and fuzziness.
        reread = SeqIO.read(base.with_suffix('.gb'), 'genbank')
        assert reread.seq == locus.seq
        assert [str(f.location) for f in reread.features] == [str(f.location) for f in locus.features]
        gff_rows = [line.split('\t') for line in base.with_suffix('.gff3').read_text(encoding='utf-8').splitlines()
                    if line and not line.startswith('#')]
        expected_gff = [('region' if f.type == 'source' else f.type, int(p.start), int(p.end),
                         {1: '+', -1: '-'}.get(p.strand, '.')) for f in locus.features for p in f.location.parts]
        assert [(r[2], int(r[3])-1, int(r[4]), r[6]) for r in gff_rows] == expected_gff
        verified.append(dict(locus_id=name, whole_RT_preserved=True, partners_preserved=True,
                            every_base_matches_source=True, fully_retained_features_checked=checked_features,
                            genbank_coordinate_roundtrip=True, gff_coordinate_roundtrip=True))
    if not records_out:
        raise ValueError('No candidate RTs to extract')
    SeqIO.write(records_out, out/'oriented_loci.fna', 'fasta')
    SeqIO.write(records_out, out/'oriented_loci.gb', 'genbank')
    write_tsv(out/'locus_manifest.tsv', manifests, list(manifests[0]))
    write_tsv(out/'coordinate_map.tsv', maps, list(maps[0]))
    write_tsv(out/'gene_annotations.tsv', genes, list(genes[0]) if genes else ['locus_id'])
    write_tsv(out/'sequence_and_annotation_issues.tsv', issues,
              ['locus_id', 'kind', 'locus_start0', 'locus_end0', 'details_json'])
    write_json(out/'RT_review.json', reviews)
    source_checks = {accession: sha((ROOT/sources[accession]['path']).read_bytes()) == sources[accession]['file_sha256']
                     for accession in cache}
    assert all(source_checks.values()), 'Original source record changed during extraction'
    write_json(out/'verification.json', dict(passed=True, loci=verified, skipped_without_RT=skipped,
                                           original_source_files_unchanged=source_checks))
    shutil.copyfile(args.candidates, out/'input_candidate_loci.tsv')
    shutil.copyfile(args.sources, out/'input_source_records.tsv')
    shutil.copyfile(__file__, out/'code'/Path(__file__).name)
    test_path = Path(__file__).with_name('test_extract_art_loci.py')
    if test_path.exists():
        shutil.copyfile(test_path, out/'code'/test_path.name)
    summary = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
        biopython=Bio.__version__, command=sys.argv, loci=len(manifests), source_records=len(cache),
        reverse_strand_RTs=sum(m['strand'] == -1 for m in manifests),
        full_upstream_flanks=sum(m['upstream_shortfall_nt'] == 0 for m in manifests),
        full_downstream_flanks=sum(m['downstream_shortfall_nt'] == 0 for m in manifests),
        wrapped_loci=sum(m['wrapped'] for m in manifests), ambiguous_nt=sum(m['ambiguous_nt'] for m in manifests),
        issues_by_kind=dict(Counter(i['kind'] for i in issues)),
        RTs_requiring_review=[r['RT_protein_id'] for r in reviews if r['review_required']],
        coordinate_convention='zero_based_half_open', gene_calling='not_used; current candidate CDSs are source annotated',
        limitations=['Extraction preserves the supplied CDS boundaries; it does not resolve missed genes or split RT annotations.',
                    'Annotation-free intervals can contain noncoding arrays; they are not inferred assembly gaps.',
                    'For unannotated contigs, run a viral gene caller and discovery before providing candidate CDSs.',
                    'Circular wrapping requires an explicit circular topology in the authoritative source record.'])
    write_json(out/'summary.json', summary)
    report = [
        '# RT-oriented locus extraction', '',
        f"Extracted {len(manifests)} candidate RT loci from {len(cache)} source records, including "
        f"{summary['reverse_strand_RTs']} reverse-strand RTs. Every locus retains its full annotated "
        f"RT and named partner CDSs. Requested flanks: {args.upstream:,} nt upstream and {args.downstream:,} nt downstream. "
        'Discovery classifications are inherited from the input candidate table and remain provisional.', '',
        '## Files', '',
        '- `oriented_loci.fna` and `oriented_loci.gb`: all oriented loci.',
        '- `loci/`: individual FASTA, GenBank, GFF3 and JSON coordinate maps.',
        '- `coordinate_map.tsv`: exact source/locus mapping segments, in locus order.',
        '- `locus_manifest.tsv`: sequence paths, RT position, partners, topology, flank availability and boundary flags.',
        '- `gene_annotations.tsv`: gene/CDS coordinates, original feature index, source location, strand and qualifiers.',
        '- `sequence_and_annotation_issues.tsv`: ambiguous sequence, explicit gaps, annotation-free intervals and clipped features.',
        '- `original_sources/`: unchanged source files and indexed annotation snapshots.',
        '- `RT_review.json`, `verification.json`, `summary.json`, `output_sha256.json`: QC and provenance.', '',
        '## Coordinate convention', '',
        'All TSV/JSON intervals are zero-based, half-open. GenBank/GFF3 are one-based, inclusive. '
        'On a mapping segment [s,e) with locus offset o, [a,b) maps to [o+a-s,o+b-s) for orientation +1, '
        'or [o+e-b,o+e-a) for orientation -1. To map a locus interval that crosses a segment boundary, '
        'split it at that boundary first. Strand in the map is the RT orientation relative to the original contig. '
        'Mapped gene strands are relative to the oriented locus.', '',
        '## Current results', '',
        f"There are {summary['ambiguous_nt']} ambiguous nucleotides and {summary['wrapped_loci']} wrapped loci. "
        f"Full requested flanks are available for {summary['full_upstream_flanks']}/{len(manifests)} upstream "
        f"and {summary['full_downstream_flanks']}/{len(manifests)} downstream intervals. "
        'Circular wrapping is enabled only by explicit circular topology in the source record; '
        'unknown/linear topology never enables wrapping.', '',
        '| Locus | Original RT strand | Length (nt) |', '|---|---:|---:|']
    report.extend(f"| {m['locus_id']} | {m['strand']:+d} | {m['length_nt']} |" for m in manifests)
    report.extend(['', f"{sum(r['source_translation_verified'] for r in reviews)}/{len(reviews)} RT translations "
        'pass checks under the annotated genetic code. Review flags:', ''])
    flagged = [r for r in reviews if r['review_required']]
    report.extend([f"- {r['RT_protein_id']}: {'; '.join(r['review_reasons'])}. "
                   f"Source translation verified: {r['source_translation_verified']}." for r in flagged] or ['- None.'])
    report.extend(['', 'Review flags retain the source CDS for inspection; profile coverage alone does not '
        'establish biological truncation. See `RT_review.json` for all sequence/annotation checks.', '',
        'All source files are unchanged. Every extracted base, fully retained feature sequence, complete RT '
        'and partner, GenBank location round trip and GFF coordinate conversion passed verification.', '',
        '## Annotation limits', '',
        'All current sources are annotated; no gene caller was used. Unannotated viral contigs require '
        'candidate CDS calls and profile discovery before this extractor can map their RTs. '
        'Pyrodigal-gv is an option: https://github.com/althonos/pyrodigal-gv . '
        'Keep caller/model/genetic-code provenance and review suspected split RTs. This step preserves '
        'the provided gene boundaries and does not repair gene calls.', '',
        'Annotation-free intervals describe annotation coverage; they may contain noncoding ART arrays '
        'and are not inferred assembly gaps. Context annotations cut by a locus boundary are marked clipped; '
        'full original qualifiers and translations remain in the archived source. Extracted clipped CDSs '
        'omit the full-length translation and adjust codon_start for the retained reading frame.', '',
        'Implementation uses Biopython 1.88 plus custom segment/interval handling. See ../../README.md '
        'and ../../scripts/extract_art_loci.py for commands and parameters.', ''])
    (out/'REPORT.md').write_text('\n'.join(report), encoding='utf-8')
    write_json(out/'output_sha256.json', {p.relative_to(out).as_posix(): sha(p.read_bytes())
        for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'output_sha256.json'})
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, default=ROOT/'results/art_discovery/candidate_loci.tsv')
    parser.add_argument('--sources', type=Path, default=ROOT/'results/art_discovery/source_records.tsv')
    parser.add_argument('--output', type=Path, default=ROOT/'results/oriented_art_loci')
    parser.add_argument('--upstream', type=int, default=6000)
    parser.add_argument('--downstream', type=int, default=6000)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
