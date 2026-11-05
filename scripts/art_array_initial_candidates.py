"""Export initial-mode candidate boundaries and flag disagreement with refinement.

No null comparisons or main-call criteria are changed. This keeps initial
detections usable for review when stringent delimitation fails or drops copies.
"""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import statistics

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from extract_art_loci import Segment, read_tsv, source_pieces, write_json, write_tsv


def overlapping_copies(first, second):
    return sum(any(max(a,x) < min(b,y) for x,y in second) for a,b in first)


def export_array_sequences(output, inputs):
    manifests = {r['locus_id']:r for r in read_tsv(inputs/'locus_manifest.tsv')}
    mapping_rows = read_tsv(inputs/'coordinate_map.tsv')
    rows = []
    for table,fasta,mode in (('arrays.tsv','arrays.fna','strict_delimitation'),
                             ('initial_arrays.tsv','initial_arrays.fna','initial_candidate')):
        records = []
        for array in read_tsv(output/table):
            locus = array['locus_id']
            record = SeqIO.read(inputs/manifests[locus]['genbank_path'],'genbank')
            a,b = int(array['array_start0']),int(array['array_end0'])
            mapping = [Segment(**{key:int(r[key]) for key in (
                'source_start0','source_end0','locus_start0','locus_end0','orientation')})
                for r in mapping_rows if r['locus_id'] == locus]
            mapping.sort(key=lambda s:s.locus_start0)
            pieces = source_pieces(a,b,mapping)
            assert sum(p['source_end0']-p['source_start0'] for p in pieces) == b-a
            records.append(SeqRecord(record.seq[a:b],id=array['array_id'],description=f'{mode}; locus=[{a},{b})'))
            rows.append(dict(array_id=array['array_id'],locus_id=locus,mode=mode,locus_start0=a,locus_end0=b,
                             length_nt=b-a,source_segments_json=json.dumps(pieces)))
        SeqIO.write(records,output/fasta,'fasta')
    write_tsv(output/'array_sequence_map.tsv',rows,list(rows[0]) if rows else ['array_id','mode'])


def export_initial_candidates(output, inputs, cfg):
    # Import lazily to keep the main scanner's process-worker import lightweight.
    from art_array_detector import align_boundaries, boundary_estimate, context_flags
    manifests = {r['locus_id']: r for r in read_tsv(inputs/'locus_manifest.tsv')}
    locus_calls = read_tsv(output/'locus_calls.tsv')
    arrays = read_tsv(output/'arrays.tsv')
    refined_repeats = read_tsv(output/'repeats.tsv')
    raw_maps = read_tsv(inputs/'coordinate_map.tsv')
    initial_arrays, repeats, spacers, review = [], [], [], []
    repeat_fasta, spacer_fasta = [], []
    extra_commands = []
    for call in locus_calls:
        locus = call['locus_id']
        flags = []
        if call['initial_pass'] != 'True':
            continue
        meta = manifests[locus]
        record = SeqIO.read(inputs/meta['genbank_path'], 'genbank')
        rt_start = int(meta['RT_locus_start0'])
        start = max(0, rt_start-cfg['initial_window_nt'])
        sequence = str(record.seq[start:rt_start]).upper()
        detail = json.loads((output/'details'/f'{locus}.scan.json').read_text())
        initial = detail['initial']
        positions = initial['positions']
        array_id = locus+'__initial_array'
        ungapped = boundary_estimate(sequence,positions,cfg['initial_seed_nt'],cfg)
        aligned = align_boundaries(sequence,positions,cfg['initial_seed_nt'],cfg,output,array_id)
        extra_commands.append(aligned['command'])
        if aligned['status'] != 'aligned' or ungapped is None:
            review.append(dict(locus_id=locus,flag='initial_boundary_estimation_failed',details=aligned['status']))
            continue
        intervals = [[a+start,b+start] for a,b in aligned['intervals']]
        context = context_flags(intervals,record.features,rt_start,cfg)
        mapping = [Segment(**{key:int(r[key]) for key in (
            'source_start0','source_end0','locus_start0','locus_end0','orientation')})
            for r in raw_maps if r['locus_id'] == locus]
        mapping.sort(key=lambda s:s.locus_start0)
        matches = []
        for array in arrays:
            if array['locus_id'] != locus:
                continue
            other = [[int(r['locus_start0']),int(r['locus_end0'])] for r in refined_repeats
                     if r['array_id'] == array['array_id']]
            if overlapping_copies(intervals,other) >= 2:
                matches.append(array)
        if not matches:
            flags.append('initial_detection_not_confirmed_by_delimitation')
        if any(int(a['copy_count']) != len(intervals) for a in matches):
            flags.append('initial_vs_refined_copy_count_disagreement')
        initial_lengths = [b-a for a,b in intervals]
        if any(statistics.median(json.loads(a['repeat_lengths_nt_json'])) != statistics.median(initial_lengths)
               for a in matches):
            flags.append('initial_vs_refined_repeat_boundary_disagreement')
        if context['coding_repeat']:
            flags.append('initial_chain_is_coding_repeat')
        if not context['adjacent_to_RT']:
            flags.append('initial_chain_not_adjacent_to_RT')
        if aligned['boundary_at_search_limit'] or ungapped['boundary_at_search_limit']:
            flags.append('initial_boundary_at_search_limit')
        starts = [a for a,b in intervals]
        lengths = [b-a for a,b in intervals]
        gaps = [c-b for (a,b),(c,d) in zip(intervals,intervals[1:])]
        initial_arrays.append(dict(array_id=array_id,locus_id=locus,accession=meta['accession'],
            discovery_label='published_style_initial_candidate', confirmation='matches_delimited_array' if matches else
            'not_confirmed_by_delimitation', copy_count=len(intervals),**context,
            repeat_consensus=aligned['consensus'], repeat_lengths_nt_json=json.dumps(lengths),
            repeat_start_spacings_nt_json=json.dumps([c-a for a,c in zip(starts,starts[1:])]),
            spacer_lengths_nt_json=json.dumps(gaps),
            matching_delimited_array_ids=';'.join(a['array_id'] for a in matches),
            matching_delimited_copy_counts=';'.join(a['copy_count'] for a in matches),
            review_flags=';'.join(flags)))
        for flag in flags:
            review.append(dict(locus_id=locus,flag=flag,details=array_id))
        for i,((a,b),aligned_text) in enumerate(zip(intervals,aligned['aligned_repeat_sequences']),1):
            repeat_id = f'{array_id}.repeat{i}'
            text = str(record.seq[a:b])
            assert text == aligned_text.replace('-','')
            segments = source_pieces(a,b,mapping)
            assert sum(p['source_end0']-p['source_start0'] for p in segments) == b-a
            repeats.append(dict(repeat_id=repeat_id,array_id=array_id,locus_id=locus,copy_index=i,
                locus_start0=a,locus_end0=b,length_nt=b-a,sequence=text,aligned_sequence=aligned_text,
                source_segments_json=json.dumps(segments)))
            repeat_fasta.append(SeqRecord(Seq(text),id=repeat_id,description=f'initial-mode candidate; locus=[{a},{b})'))
            if i < len(intervals):
                c = intervals[i][0]
                spacer_id = f'{array_id}.spacer{i}'
                text = str(record.seq[b:c])
                assert c-b == (c-a)-(b-a)
                spacers.append(dict(spacer_id=spacer_id,array_id=array_id,locus_id=locus,spacer_index=i,
                    locus_start0=b,locus_end0=c,length_nt=c-b,repeat_start_spacing_nt=c-a,
                    preceding_repeat_length_nt=b-a,sequence=text,
                    source_segments_json=json.dumps(source_pieces(b,c,mapping))))
                spacer_fasta.append(SeqRecord(Seq(text),id=spacer_id,description=f'initial-mode candidate; locus=[{b},{c})'))
    for name,rows,empty in (('initial_arrays.tsv',initial_arrays,['array_id','locus_id']),
        ('initial_repeats.tsv',repeats,['repeat_id','array_id']),
        ('initial_spacers.tsv',spacers,['spacer_id','array_id']),
        ('review_flags.tsv',review,['locus_id','flag','details'])):
        write_tsv(output/name,rows,list(rows[0]) if rows else empty)
    SeqIO.write(repeat_fasta,output/'initial_repeats.fna','fasta')
    SeqIO.write(spacer_fasta,output/'initial_spacers.fna','fasta')
    export_array_sequences(output,inputs)
    commands_path = output/'commands.json'
    commands = json.loads(commands_path.read_text()) if commands_path.exists() else []
    initial_ids = {str(output/'alignments'/f"{a['array_id']}.mafft.fna") for a in initial_arrays}
    commands = [c for c in commands if not any(p in c for p in initial_ids)]
    write_json(commands_path,commands+extra_commands)
    summary = json.loads((output/'summary.json').read_text())
    summary.update(initial_candidate_arrays=len(initial_arrays),initial_candidate_repeat_copies=len(repeats),
                   initial_candidates_not_confirmed=sum(not a['matching_delimited_array_ids'] for a in initial_arrays),
                   loci_with_mode_disagreement=len({r['locus_id'] for r in review if 'disagreement' in r['flag']}))
    write_json(output/'summary.json',summary)
    verification_path = output/'verification.json'
    verification = json.loads(verification_path.read_text())
    verification.update(initial_candidate_sequence_coordinates_checked=True,
                        initial_candidate_spacer_formula_checked=True)
    write_json(verification_path,verification)
    shutil.copyfile(__file__,output/'code'/Path(__file__).name)
    with (output/'REPORT.md').open('a',encoding='utf-8') as handle:
        handle.write('\n## Initial candidates and disagreement review\n\n'
            '`initial_arrays.tsv`, `initial_repeats.tsv`, `initial_spacers.tsv` and matching FASTAs '
            'retain initial-mode candidates even when extended delimitation fails. Their boundaries '
            'are separately estimated by the same ungapped/MAFFT method. These are candidate annotations, '
            'and are not promoted to confirmed delimitation calls. `review_flags.tsv` identifies '
            'unconfirmed initial detections and copy-count/boundary disagreements.\n\n'
            '`arrays.fna` and `initial_arrays.fna` contain complete candidate array tracts, '
            'including their spacers. `array_sequence_map.tsv` maps these full tracts to source coordinates.\n\n'
            '| Locus | Initial copies | Refined copies at same tract | Review flags |\n'
            '|---|---:|---|---|\n')
        for row in initial_arrays:
            handle.write(f"| {row['locus_id']} | {row['copy_count']} | {row['matching_delimited_copy_counts'] or 'none'} | "
                         f"{row['review_flags'] or 'none'} |\n")
        handle.write('\nRepeat boundaries and copy counts remain computational estimates. A strict PWM can '
            'retain a more conserved subset of an array and expand the inferred boundary into flanks shared '
            'by that subset. Disagreement with the initial seed chain requires manual review before '
            'using spacer sequences for downstream structural or evolutionary inference.\n')
    return summary
