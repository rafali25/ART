"""ART array inference: fresh-seed nulls, dinucleotides, cohort correction and QC.

Consumes a checksum-verified replication scan and retains its paper-style labels.
All coordinates are zero-based, half-open. No universal repeat sequence is used.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import sys
import time

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from art_array_detector import (ROOT, BASES, words, decode, approximate_positions,
    best_chain, boundary_estimate, base_background, refinement_scan, shuffled,
    stable_seed, empirical, pwm_scores, nonoverlap_scored, align_boundaries,
    context_flags, initial_support, distinct_arrays, conserved_block, dna_array)
from art_array_initial_candidates import overlapping_copies
from extract_art_loci import read_tsv, write_tsv, write_json, Segment, source_pieces


@lru_cache(maxsize=128)
def dinucleotide_graph(sequence):
    """Enumerate edge-weighted last-exit trees, rooted at the final base.

    There are at most 4 DNA vertices and 4^3 candidate parent assignments.
    Tree weight is the product of edge multiplicities. This samples labelled
    arborescences uniformly, then last-exit Euler trails uniformly (the Euler
    algorithm). Enumeration replaces Wilson sampling for this tiny alphabet.
    """
    vertices = sorted(set(sequence))
    root = sequence[-1]
    edges = {v:[] for v in vertices}
    for a,b in zip(sequence,sequence[1:]):
        edges[a].append(b)
    nonroot = [v for v in vertices if v != root]
    counts = {v:Counter(edges[v]) for v in vertices}
    choices = [[v for v in vertices if counts[u][v]] for u in nonroot]
    trees,weights = [],[]
    for parents in itertools.product(*choices):
        tree = dict(zip(nonroot,parents))
        valid = True
        for v in nonroot:
            seen = set()
            while v != root:
                if v in seen or v not in tree:
                    valid = False
                    break
                seen.add(v)
                v = tree[v]
            if not valid:
                break
        if valid:
            trees.append(tree)
            weights.append(math.prod(counts[v][tree[v]] for v in nonroot))
    if not trees:
        raise ValueError('No endpoint-rooted last-exit tree for valid sequence')
    probabilities = np.array(weights,dtype=float)
    probabilities /= probabilities.sum()
    return edges,trees,probabilities


def dinucleotide_shuffle(sequence, rng):
    """Exact mono/di counts and endpoints, with ambiguity positions unchanged.

    ACGT runs are randomized separately; fixed endpoints also preserve pairs
    across ambiguous symbols. This intentionally preserves more than global
    dinucleotide composition when ambiguity separates runs.
    """
    import re
    result = list(sequence)
    for match in re.finditer('[ACGT]+',sequence):
        text = match.group()
        if len(text) < 3:
            continue
        edges,trees,probabilities = dinucleotide_graph(text)
        tree = trees[int(rng.choice(len(trees),p=probabilities))]
        lists = {}
        for v,neighbours in edges.items():
            remaining = list(neighbours)
            last = tree.get(v)
            if last is not None:
                remaining.remove(last)
            rng.shuffle(remaining)
            if last is not None:
                remaining.append(last)
            lists[v] = remaining
        cursor = {v:0 for v in edges}
        output = [text[0]]
        for i in range(len(text)-1):
            v = output[-1]
            output.append(lists[v][cursor[v]])
            cursor[v] += 1
        assert output[-1] == text[-1]
        assert all(cursor[v] == len(edges[v]) for v in edges)
        result[match.start():match.end()] = output
    return ''.join(result)


def sensitivity_batch(task):
    sequence,kind,cfg,seeds = task
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        permuted = dinucleotide_shuffle(sequence,rng) if kind == 'dinucleotide' else shuffled(sequence,rng,cfg['shuffle_block_nt'])
        # Full rediscovery and all-seed best-chain maximization, never a fixed PWM.
        result = refinement_scan(permuted,cfg,False)
        rows.append((seed,result['score']))
    return rows


def inferential_null(sequence,kind,cfg,assessment,n,locus,window,pool,existing=()):
    rows = list(existing[:n])
    start = len(rows)
    if kind == 'local_blocks':
        seeds = [stable_seed(cfg['random_seed'],locus,'refinement_'+window,len(sequence),i) for i in range(start,n)]
    else:
        seeds = [stable_seed(assessment['random_seed'],locus,kind,window,i) for i in range(start,n)]
    jobs = [(sequence,kind,cfg,seeds[i:i+20]) for i in range(0,len(seeds),20)]
    for batch in pool.map(sensitivity_batch,jobs):
        rows.extend(batch)
    return rows


def bonferroni(p, loci):
    if loci < 1:
        raise ValueError('The testing family must include at least one locus')
    return min(1.,float(p)*loci)


def entropy(sequence):
    counts = Counter(b for b in sequence if b in BASES)
    n = sum(counts.values())
    return -sum((count/n)*math.log2(count/n) for count in counts.values()) if n else 0.


def short_period_scan(sequence,cfg,max_period=59):
    k = cfg['refinement_seed_nt']
    index = words(sequence,k)
    best = dict(score=0.,copy_count=0,positions=[],base_spacing=None,seed='')
    for code,exact in index.items():
        if len(exact) < 2:
            continue
        positions = approximate_positions(index,code,k,cfg['refinement_max_mismatches'])
        chain = best_chain(positions,k,max_period,cfg['refinement_spacing_tolerance'],strict_gaps=True)
        if len(chain['positions']) < 3:
            continue
        boundary = boundary_estimate(sequence,chain['positions'],k,cfg)
        if boundary is None:
            continue
        score = (len(chain['positions'])-1)*boundary['information_bits']
        if score > best['score']:
            best = dict(score=score,copy_count=len(chain['positions']),positions=chain['positions'],
                        base_spacing=chain['base_spacing'],seed=decode(code,k))
    return best


def quality_flags(record,intervals,rt_start,meta,cfg,assessment,short,chain_score,initial_copies=None):
    import re
    context = context_flags(intervals,record.features,rt_start,cfg)
    a,b = intervals[0][0],intervals[-1][1]
    repeat_texts = [str(record.seq[x:y]).upper() for x,y in intervals]
    array = str(record.seq[a:b]).upper()
    overlaps = set()
    intervening = set()
    cds_names = {f.qualifiers['locus_tag'][0]:f.qualifiers.get('protein_id',f.qualifiers['locus_tag'])[0]
                 for f in record.features if f.type == 'CDS' and f.qualifiers.get('locus_tag')}
    for f in record.features:
        if f.type not in ('CDS','gene','tRNA','rRNA','ncRNA','tmRNA') or f.location is None:
            continue
        identifier = f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',f.qualifiers.get('gene',['unknown'])))[0]
        identifier = cds_names.get(identifier,identifier)
        if f.type == 'CDS' and any(max(x,int(p.start)) < min(y,int(p.end)) for x,y in intervals for p in f.location.parts):
            overlaps.add(identifier)
        if b <= int(f.location.start) and int(f.location.end) <= rt_start:
            intervening.add(identifier)
    low = any(entropy(text) < assessment['low_complexity_entropy_bits'] or
              re.search(r'(.)\1{'+str(assessment['low_complexity_homopolymer_nt']-1)+',}',text)
              for text in repeat_texts)
    flags = []
    if overlaps: flags.append('repeat_overlaps_CDS')
    if context['coding_repeat']: flags.append('triplet_period_coding_repeat')
    if intervening: flags.append('gene_between_array_and_RT')
    if low: flags.append('low_complexity_repeat')
    if entropy(array) < assessment['low_complexity_entropy_bits']: flags.append('low_complexity_array_region')
    ext_start = max(0,rt_start-cfg['extended_window_nt'])
    if short['score'] > chain_score:
        if sum(a <= p+ext_start < b for p in short['positions']) >= 3:
            flags.append('stronger_short_period_pattern')
        else:
            flags.append('stronger_short_period_elsewhere_in_locus')
    if any(c not in BASES for c in array): flags.append('ambiguous_bases_in_array')
    upstream_ambiguity = sum(c not in BASES for c in str(record.seq[ext_start:rt_start]).upper())
    if upstream_ambiguity: flags.append('ambiguous_bases_in_search_window')
    if int(meta.get('upstream_shortfall_nt',0)) > 0: flags.append('contig_limited_upstream')
    if meta.get('source_origin_reached') == 'True' or meta.get('source_end_reached') == 'True':
        flags.append('locus_reaches_source_boundary')
    if a <= ext_start: flags.append('array_touches_upstream_search_boundary')
    if initial_copies is not None and initial_copies != len(intervals): flags.append('initial_vs_final_copy_count_disagreement')
    return dict(flags=flags,context=context,CDS_overlap_ids=sorted(overlaps),
                intervening_gene_ids=sorted(intervening),repeat_entropies=[entropy(t) for t in repeat_texts],
                array_entropy=entropy(array),ambiguous_array_nt=sum(c not in BASES for c in array),
                ambiguous_upstream_nt=upstream_ambiguity)


def main_status(copies,adjusted_p,flags,alpha=.05):
    if copies == 2:
        return 'exploratory_two_copy_not_accepted'
    if copies < 3:
        return 'insufficient_recovered_copies'
    if adjusted_p > alpha:
        return 'not_supported_after_locus_correction'
    blocking = {'repeat_overlaps_CDS','triplet_period_coding_repeat','low_complexity_repeat','low_complexity_array_region',
        'stronger_short_period_pattern','ambiguous_bases_in_array','array_touches_upstream_search_boundary',
        'gene_between_array_and_RT','initial_vs_final_copy_count_disagreement','boundary_search_limit'}
    if blocking.intersection(flags):
        return 'statistically_supported_requires_review'
    return 'accepted_array_candidate'


def aligned_pwm(sequence,chain,cfg,output,array_id):
    """Align first, estimate boundaries, THEN build PWM and search for copies."""
    first = align_boundaries(sequence,chain['positions'],cfg['refinement_seed_nt'],cfg,output,array_id+'.training')
    if first['status'] != 'aligned':
        return dict(status=first['status'],training=first)
    texts = first['aligned_repeat_sequences']
    # Explicitly project onto columns occupied in >=80% of aligned copies. Gaps
    # in minority copies are missing observations, not an extra DNA character.
    # The original gapped alignment and projection mask remain in the evidence.
    columns = [j for j in range(len(texts[0])) if sum(t[j] in BASES for t in texts) >=
               cfg['boundary_consensus_fraction']*len(texts)-1e-9]
    if not columns:
        return dict(status='no_occupied_alignment_columns',training=first)
    counts = np.array([[sum(t[j] == b for t in texts) for b in BASES] for j in columns],dtype=float)
    probabilities = (counts+cfg['pwm_pseudocount_per_base'])/(counts.sum(axis=1,keepdims=True)+4*cfg['pwm_pseudocount_per_base'])
    pwm = np.log2(probabilities/base_background(sequence))
    null = []
    for i in range(cfg['pwm_shuffles']):
        seed = stable_seed(cfg['random_seed'],array_id,'aligned_PWM',i)
        null.append((seed,float(pwm_scores(shuffled(sequence,np.random.default_rng(seed),cfg['shuffle_block_nt']),pwm).max())))
    threshold = max(score for _,score in null)
    scores = pwm_scores(sequence,pwm)
    hits = nonoverlap_scored(np.flatnonzero(scores > threshold),scores,len(pwm))
    selected = best_chain(hits,*cfg['refinement_base_spacing_nt'],cfg['refinement_spacing_tolerance'],
                          missing=cfg['maximum_missing_copies_per_chain'] if cfg['allow_single_missing_copy'] else 0)
    intervals = [[p,p+len(pwm)] for p in selected['positions']]
    final_sequences = [sequence[a:b] for a,b in intervals]
    if len(intervals) >= 2:
        matrix = np.array([dna_array(t) for t in final_sequences])
        block = conserved_block(matrix,len(pwm)//2,cfg['boundary_consensus_fraction'],cfg['boundary_allowed_lapses'])
        if block is None:
            return dict(status='no_conserved_block_after_PWM',training=first,pwm=pwm.tolist())
        a,b = block['start'],block['end']
        intervals = [[p+a,p+b] for p in selected['positions']]
        final_sequences = [sequence[a:b] for a,b in intervals]
    return dict(status='rescanned',training=first,retained_alignment_columns=columns,
        alignment_projection_used=len(columns) != len(texts[0]),
        training_alignment_has_indels=any('-' in t for t in texts),
        pwm=pwm.tolist(),PWM_null=null,PWM_threshold=threshold,
        all_significant_starts=hits,selected=selected,intervals=[[p,p+len(pwm)] for p in selected['positions']],
        final_intervals=intervals,final_aligned_sequences=final_sequences,
        PWM_scores=[float(scores[p]) for p in selected['positions']])


def _scan_core_hash():
    import inspect
    functions = (refinement_scan,words,decode,approximate_positions,best_chain,boundary_estimate,
                 base_background,conserved_block,dna_array,shuffled,dinucleotide_graph.__wrapped__,dinucleotide_shuffle)
    return hashlib.sha256('\n'.join(inspect.getsource(f) for f in functions).encode()).hexdigest()


def run(args):
    started = time.perf_counter()
    baseline,inputs,output = args.baseline.resolve(),args.loci.resolve(),args.output.resolve()
    cfg = json.loads((baseline/'run_config.json').read_text())['art_array_search']
    assessment = json.loads(args.config.read_text())
    if cfg['minimum_copies_for_main_call'] < 3 or assessment['minimum_main_copies'] < 3:
        raise ValueError('Main evidence requires >=3 copies; two-copy candidates are exploratory only')
    if assessment['multiple_testing_method'] != 'bonferroni_across_all_input_loci':
        raise ValueError('Unsupported multiple-testing method')
    for folder in (output,output/'nulls',output/'alignments',output/'details',output/'code'):
        folder.mkdir(parents=True,exist_ok=True)
    for path,digest in json.loads((baseline/'output_sha256.json').read_text()).items():
        if hashlib.sha256((baseline/path).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Baseline artifact changed: {path}')
    provenance = json.loads((baseline/'provenance.json').read_text())
    for path,digest in provenance['input_hashes'].items():
        if hashlib.sha256((inputs/path).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Baseline locus changed: {path}')
    manifests = read_tsv(inputs/'locus_manifest.tsv')
    families = len(manifests)
    if not families:
        raise ValueError('No loci in testing family')
    old_nulls = defaultdict(list)
    for row in read_tsv(baseline/'null_scores.tsv'):
        if row['stage'].startswith('refinement_local_blocks_'):
            old_nulls[row['locus_id'],row['stage'].rsplit('_',1)[-1]].append((int(row['rng_seed']),float(row['score'])))
    locus_rows,null_rows = [],[]
    null_by_locus = {}
    core_hash = _scan_core_hash()
    with ProcessPoolExecutor(max_workers=assessment['workers']) as pool:
        for i,meta in enumerate(manifests,1):
            locus = meta['locus_id']
            record = SeqIO.read(inputs/meta['genbank_path'],'genbank')
            rt = int(meta['RT_locus_start0'])
            sequences = {'3000':str(record.seq[max(0,rt-cfg['initial_window_nt']):rt]).upper(),
                         '6000':str(record.seq[max(0,rt-cfg['extended_window_nt']):rt]).upper()}
            scan = json.loads((baseline/'details'/f'{locus}.scan.json').read_text())
            cache_path = output/'nulls'/f'{locus}.json'
            key = dict(core_sha256=core_hash,cfg=cfg,assessment=assessment,
                       sequence_sha256={k:hashlib.sha256(s.encode()).hexdigest() for k,s in sequences.items()})
            if args.reuse_nulls and cache_path.exists():
                cached = json.loads(cache_path.read_text())
                if cached['key'] != key:
                    raise ValueError(f'Cannot reuse incompatible null cache: {locus}')
                rows_by_test = cached['scores']
                print(f'[{i}/{families}] {locus}: verified cached nulls',flush=True)
            else:
                rows_by_test = {}
                for kind,count in (('local_blocks',assessment['local_inference_shuffles_per_window']),
                                   ('dinucleotide',assessment['dinucleotide_shuffles_per_window'])):
                    for window,sequence in sequences.items():
                        existing = old_nulls[locus,window] if kind == 'local_blocks' else []
                        print(f'[{i}/{families}] {locus}: {kind}, {window} nt, {count} fixed draws',flush=True)
                        rows_by_test[kind+'_'+window] = inferential_null(sequence,kind,cfg,assessment,count,locus,window,pool,existing)
                write_json(cache_path,dict(key=key,scores=rows_by_test))
            null_by_locus[locus] = rows_by_test
            pvalues = {}
            row = dict(locus_id=locus,testing_family_loci=families,
                       inferential_window='6000',diagnostic_window='3000',
                       actual_initial_window_nt=len(sequences['3000']),actual_extended_window_nt=len(sequences['6000']))
            for kind in ('local_blocks','dinucleotide'):
                for window in ('3000','6000'):
                    name = kind+'_'+window
                    observed = scan['refinement_'+window]['score']
                    test = empirical(observed,[score for seed,score in rows_by_test[name]])
                    pvalues[name] = test['empirical_p']
                    row[name+'_best_score'] = observed
                    row[name+'_p'] = test['empirical_p']
                    row[name+'_null_max'] = test['null_max']
                    row[name+'_shuffles'] = test['shuffles']
                    for j,(seed,score) in enumerate(rows_by_test[name]):
                        null_rows.append(dict(locus_id=locus,background=kind,window_nt=window,
                            shuffle_index=j,rng_seed=seed,best_chain_score=score))
            inferential_p = max(p for name,p in pvalues.items() if name.endswith('_6000'))
            row.update(max_background_p=inferential_p,
                       bonferroni_p=bonferroni(inferential_p,families),
                       statistic='best_all_seed_chain_score_per_window',
                       whole_seed_discovery_rerun_in_every_null=True)
            locus_rows.append(row)
    write_tsv(output/'background_tests.tsv',locus_rows,list(locus_rows[0]))
    write_tsv(output/'null_best_scores.tsv',null_rows,list(null_rows[0]))
    produce_candidates(output,inputs,baseline,cfg,assessment,manifests,null_by_locus,locus_rows)
    summary = json.loads((output/'summary.json').read_text())
    summary.update(elapsed_seconds=round(time.perf_counter()-started,2),
        fixed_null_draws_per_window=dict(local_blocks=assessment['local_inference_shuffles_per_window'],
                                       dinucleotide=assessment['dinucleotide_shuffles_per_window']),
        calibration_draws_including_reused=len(null_rows),testing_family_loci=families,
        multiple_testing_method=assessment['multiple_testing_method'],familywise_alpha=assessment['familywise_alpha'])
    write_json(output/'summary.json',summary)
    write_json(output/'run_config.json',dict(search=cfg,assessment=assessment))
    write_json(output/'provenance.json',dict(command=sys.argv,baseline=str(baseline),input_loci=str(inputs),
        core_sha256=core_hash,baseline_checksum_manifest_sha256=hashlib.sha256((baseline/'output_sha256.json').read_bytes()).hexdigest(),
        reused_local_nulls='fixed-length prefix of verified baseline draws; missing draws generated with identical seed schedule',
        dinucleotide_method_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/',
        dinucleotide_algorithm='Euler last-exit trees; exact weighted enumeration for <=4 DNA vertices'))
    for name in ('art_array_assessment.py','art_array_detector.py','art_array_initial_candidates.py',
                 'extract_art_loci.py','test_art_array_assessment.py'):
        path = Path(__file__).with_name(name)
        if path.exists(): shutil.copyfile(path,output/'code'/name)
    report(output,summary,locus_rows,assessment)
    write_json(output/'output_sha256.json',{p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'output_sha256.json'})
    print(json.dumps(summary,indent=2),flush=True)


def produce_candidates(output,inputs,baseline,cfg,assessment,manifests,null_by_locus,locus_rows):
    maps_raw = read_tsv(inputs/'coordinate_map.tsv')
    initial = {r['locus_id']:r for r in read_tsv(baseline/'initial_arrays.tsv')}
    old_strict = {r['locus_id']:r for r in read_tsv(baseline/'arrays.tsv')}
    array_rows,repeat_rows,spacer_rows,two_rows = [],[],[],[]
    repeat_seqs,spacer_seqs,array_seqs = [],[],[]
    for meta in manifests:
        locus = meta['locus_id']
        record = SeqIO.read(inputs/meta['genbank_path'],'genbank')
        rt = int(meta['RT_locus_start0'])
        offset = max(0,rt-cfg['extended_window_nt'])
        sequence = str(record.seq[offset:rt]).upper()
        scan = json.loads((baseline/'details'/f'{locus}.scan.json').read_text())
        candidates = scan['refinement_6000']['candidates']
        short = short_period_scan(sequence,cfg,assessment['short_period_max_nt'])
        mapping = [Segment(**{k:int(r[k]) for k in ('source_start0','source_end0','locus_start0','locus_end0','orientation')})
                   for r in maps_raw if r['locus_id'] == locus]
        mapping.sort(key=lambda s:s.locus_start0)
        # All seed chains retain their original scan records. Expensive alignment
        # is limited to candidates whose full-search p can pass the cohort rule.
        eligible = []
        for candidate in candidates:
            p = max(empirical(candidate['score'],[s for _,s in values])['empirical_p']
                    for name,values in null_by_locus[locus].items() if name.endswith('_6000'))
            candidate['_p'] = p
            if bonferroni(p,len(manifests)) <= assessment['familywise_alpha']:
                eligible.append(candidate)
        # Retain the best scanned candidate even when it fails, for assessment.
        if not eligible and candidates:
            eligible = [candidates[0]]
        for j,candidate in enumerate(distinct_arrays(eligible),1):
            array_id = f'{locus}__assessed{j}'
            print(f'  {array_id}: align repeat-centred segments, fit PWM, rescan',flush=True)
            model = aligned_pwm(sequence,candidate,cfg,output,array_id)
            actual_intervals = model.get('final_intervals',[])
            recovered_copies = len(actual_intervals)
            raw_intervals = actual_intervals or candidate['boundary']['intervals']
            intervals = [[a+offset,b+offset] for a,b in raw_intervals]
            tests = {name:empirical(candidate['score'],[s for _,s in values]) for name,values in null_by_locus[locus].items()}
            p = max(t['empirical_p'] for name,t in tests.items() if name.endswith('_6000'))
            adjusted = bonferroni(p,len(manifests))
            old_initial = initial.get(locus)
            initial_copies = None
            if old_initial:
                old_repeats = [[int(r['locus_start0']),int(r['locus_end0'])] for r in read_tsv(baseline/'initial_repeats.tsv')
                               if r['locus_id'] == locus]
                if overlapping_copies(intervals,old_repeats) >= 2:
                    initial_copies = int(old_initial['copy_count'])
            quality = quality_flags(record,intervals,rt,meta,cfg,assessment,short,candidate['score'],initial_copies)
            flags = quality['flags']
            if model['status'] != 'rescanned': flags.append(model['status'])
            if candidate['boundary']['boundary_at_search_limit'] or model.get('training',{}).get('boundary_at_search_limit'):
                flags.append('boundary_search_limit')
            if model['status'] == 'rescanned' and recovered_copies == 2:
                status = 'exploratory_two_copy_not_accepted'
            elif adjusted > assessment['familywise_alpha']:
                status = 'not_supported_after_locus_correction'
            elif model['status'] == 'rescanned' and recovered_copies < assessment['minimum_main_copies']:
                status = 'insufficient_recovered_copies'
            elif model['status'] != 'rescanned':
                status = 'alignment_or_PWM_requires_review'
            else:
                status = main_status(recovered_copies,adjusted,flags,assessment['familywise_alpha'])
            lengths = [b-a for a,b in intervals]
            spacings = [c-a for (a,b),(c,d) in zip(intervals,intervals[1:])]
            spacers = [c-b for (a,b),(c,d) in zip(intervals,intervals[1:])]
            info = 0.
            if model.get('final_aligned_sequences'):
                texts = model['final_aligned_sequences']
                counts = np.array([[sum(t[k] == base for t in texts) for base in BASES] for k in range(len(texts[0]))],dtype=float)
                freq = counts/counts.sum(axis=1,keepdims=True)
                info = float(np.sum(np.where(freq > 0,freq*np.log2(np.maximum(freq,1e-15)/base_background(sequence)),0)))
            row = dict(array_id=array_id,locus_id=locus,accession=meta['accession'],RT_protein_id=meta['RT_protein_id'],
                assessment_status=status,copy_count=recovered_copies,
                coordinate_stage='PWM_recovered' if actual_intervals else 'undelimited_seed_chain',
                seed_chain_copies=len(candidate['positions']),seed=candidate['seed'],seed_chain_score=candidate['score'],
                complete_array_score_descriptive=(len(intervals)-1)*info,
                statistical_score_type='seed_chain_score_tested_against_rediscovered_null_maxima',
                max_background_p=p,bonferroni_p=adjusted,testing_family_loci=len(manifests),
                paper_style_replication_delimited=locus in old_strict,
                initial_mode_pass=scan['initial']['R'] >= cfg['minimum_copies_for_main_call'] and
                                  scan['initial_test']['exceeds_all_shuffles'] and not scan['initial']['dense_veto'],
                initial_copy_count=initial_copies,flags=';'.join(flags),
                underlying_spacing_nt=model.get('selected',{}).get('base_spacing'),
                inferred_skipped_copy_links=model.get('selected',{}).get('missing_copies'),
                coding_overlap_ids=';'.join(quality['CDS_overlap_ids']),intervening_gene_ids=';'.join(quality['intervening_gene_ids']),
                stronger_short_period_score=short['score'],repeat_entropies_json=json.dumps(quality['repeat_entropies']),
                ambiguous_array_nt=quality['ambiguous_array_nt'],ambiguous_upstream_nt=quality['ambiguous_upstream_nt'],
                array_start0=intervals[0][0],array_end0=intervals[-1][1],array_span_nt=intervals[-1][1]-intervals[0][0],
                source_segments_json=json.dumps(source_pieces(intervals[0][0],intervals[-1][1],mapping)),
                seed_intervals0_json=json.dumps([[a+offset,b+offset] for a,b in candidate['boundary']['intervals']]),
                repeat_lengths_nt_json=json.dumps(lengths if actual_intervals else []),
                repeat_start_spacings_nt_json=json.dumps(spacings if actual_intervals else []),
                spacer_lengths_nt_json=json.dumps(spacers if actual_intervals else []))
            for name,test in tests.items():
                row[name+('_p' if name.endswith('_6000') else '_diagnostic_tail_fraction')] = test['empirical_p']
            array_rows.append(row)
            write_json(output/'details'/f'{array_id}.json',dict(summary=row,aligned_PWM=model,chain=candidate,
                background_tests=tests,short_period=short,quality=quality))
            if model['status'] == 'rescanned' and actual_intervals:
                aligned_records = []
                for i,(a,b) in enumerate(intervals,1):
                    repeat_id = f'{array_id}.repeat{i}'
                    text = str(record.seq[a:b])
                    original = source_pieces(a,b,mapping)
                    assert sum(p['source_end0']-p['source_start0'] for p in original) == b-a
                    repeat_rows.append(dict(repeat_id=repeat_id,array_id=array_id,locus_id=locus,copy_index=i,
                        locus_start0=a,locus_end0=b,length_nt=b-a,sequence=text,source_segments_json=json.dumps(original)))
                    item = SeqRecord(Seq(text),id=repeat_id,description=status)
                    repeat_seqs.append(item)
                    aligned_records.append(item)
                    if i < len(intervals):
                        c = intervals[i][0]
                        text = str(record.seq[b:c])
                        spacer_id = f'{array_id}.spacer{i}'
                        assert c-b == (c-a)-(b-a)
                        spacer_rows.append(dict(spacer_id=spacer_id,array_id=array_id,locus_id=locus,spacer_index=i,
                            locus_start0=b,locus_end0=c,length_nt=c-b,sequence=text,
                            repeat_start_spacing_nt=c-a,preceding_repeat_length_nt=b-a,
                            possible_missing_repeat=(c-a > model['selected']['base_spacing']*(1+cfg['refinement_spacing_tolerance'])+1e-8),
                            source_segments_json=json.dumps(source_pieces(b,c,mapping))))
                        spacer_seqs.append(SeqRecord(Seq(text),id=spacer_id,description=status))
                SeqIO.write(aligned_records,output/'alignments'/f'{array_id}.final_repeat_alignment.fna','fasta')
                array_seqs.append(SeqRecord(record.seq[intervals[0][0]:intervals[-1][1]],id=array_id,description=status))
        # Independent exploratory scan; it never contributes an accepted call or
        # borrows >=3-copy p-values. All recurring seeds are still considered.
        two_cfg = dict(cfg,minimum_copies_for_main_call=2)
        twos = refinement_scan(sequence,two_cfg)['candidates']
        pairs = distinct_arrays([c for c in twos if len(c['positions']) == 2])
        pairs.sort(key=lambda c:-c['score'])
        for j,candidate in enumerate(pairs[:assessment['exploratory_two_copy_per_locus']],1):
            array_id = f'{locus}__two_copy{j}'
            intervals = [[a+offset,b+offset] for a,b in candidate['boundary']['intervals']]
            quality = quality_flags(record,intervals,rt,meta,cfg,assessment,short,candidate['score'])
            gap = intervals[1][0]-intervals[0][1]
            if gap < 0: quality['flags'].append('overlapping_estimated_repeat_boundaries')
            two_rows.append(dict(array_id=array_id,locus_id=locus,copy_count=2,
                assessment_status='exploratory_two_copy_not_accepted',score=candidate['score'],
                calibrated_p='',bonferroni_p='',flags=';'.join(quality['flags']),
                repeat_intervals0_json=json.dumps(intervals),
                source_segments_json=json.dumps([source_pieces(a,b,mapping) for a,b in intervals]),
                repeat_sequences_json=json.dumps([str(record.seq[a:b]) for a,b in intervals]),
                spacer_sequence=str(record.seq[intervals[0][1]:intervals[1][0]]) if gap >= 0 else '',
                spacer_length_nt=gap if gap >= 0 else '',
                spacer_status='extracted' if gap >= 0 else 'undefined_overlapping_repeat_bounds',
                boundary_method='ungapped_seed_centred_exploratory_not_well_delimited'))
    for filename,rows,empty in (('array_assessments.tsv',array_rows,['array_id','assessment_status']),
        ('repeats.tsv',repeat_rows,['repeat_id','array_id']),('spacers.tsv',spacer_rows,['spacer_id','array_id']),
        ('exploratory_two_copy.tsv',two_rows,['array_id','assessment_status'])):
        write_tsv(output/filename,rows,list(rows[0]) if rows else empty)
    SeqIO.write(repeat_seqs,output/'repeats.fna','fasta')
    SeqIO.write(spacer_seqs,output/'spacers.fna','fasta')
    SeqIO.write(array_seqs,output/'arrays.fna','fasta')
    locus_assessments = []
    for meta in manifests:
        locus = meta['locus_id']
        statuses = [r['assessment_status'] for r in array_rows if r['locus_id'] == locus]
        if 'accepted_array_candidate' in statuses:
            status = 'accepted_array_candidate'
        elif 'statistically_supported_requires_review' in statuses:
            status = 'statistically_supported_requires_review'
        elif statuses:
            status = 'no_accepted_array_candidate'
        else:
            status = 'no_eligible_array_candidate'
        background = next(r for r in locus_rows if r['locus_id'] == locus)
        locus_assessments.append(dict(locus_id=locus,assessment_status=status,
            accepted_arrays=statuses.count('accepted_array_candidate'),
            assessed_arrays=len(statuses),bonferroni_best_chain_p=background['bonferroni_p'],
            exploratory_two_copy_candidates=sum(r['locus_id'] == locus for r in two_rows)))
    write_tsv(output/'locus_assessments.tsv',locus_assessments,list(locus_assessments[0]))
    accepted_ids = {r['array_id'] for r in array_rows if r['assessment_status'] == 'accepted_array_candidate'}
    SeqIO.write([r for r in array_seqs if r.id in accepted_ids],output/'accepted_arrays.fna','fasta')
    SeqIO.write([r for r in repeat_seqs if r.id.rsplit('.repeat',1)[0] in accepted_ids],output/'accepted_repeats.fna','fasta')
    SeqIO.write([r for r in spacer_seqs if r.id.rsplit('.spacer',1)[0] in accepted_ids],output/'accepted_spacers.fna','fasta')
    summary = dict(loci_scanned=len(manifests),assessed_arrays=len(array_rows),accepted_arrays=len(accepted_ids),
        assessment_counts=dict(Counter(r['assessment_status'] for r in array_rows)),
        exploratory_two_copy_candidates=len(two_rows),exported_repeat_copies=len(repeat_rows),exported_spacers=len(spacer_rows),
        coordinates='zero_based_half_open',no_fixed_MarsHill_motif=True,
        alignment_precedes_PWM=True,full_seed_search_in_all_local_and_dinucleotide_nulls=True)
    write_json(output/'summary.json',summary)


def report(output,summary,loci,assessment):
    arrays = read_tsv(output/'array_assessments.tsv')
    lines = ['# Tolerant ART array discovery with cohort assessment','',
        'The published-style replication results remain unchanged in `../art_arrays/`. This assessment '
        'adds a fixed-size inferential local-block null, an exact dinucleotide null and a correction across '
        'every input locus, including those without candidates. No fixed repeat motif is used.','',
        '## Acceptance rule','',
        f"Main candidates require >=3 PWM-recovered copies and Bonferroni-adjusted p <= {assessment['familywise_alpha']}. "
        'Candidate p is the maximum of empirical tail probabilities against the two full 6-kb-search nulls: '
        'local-block and dinucleotide. Each uses '
        '(1 + null best scores >= observed chain score)/(1 + fixed draws). Multiplying the maximum p '
        'by all input loci controls the family across loci without an independence assumption. '
        'The 6-kb null maximizes over all rediscovered seeds/chains, so candidate scores are compared '
        'with a within-window search maximum rather than a fitted-motif null. The separate per-locus '
        '`background_tests.tsv` also tests the 3-kb window using its own observed best score as a '
        'diagnostic. Candidate 6-kb scores compared with 3-kb shuffle maxima are labelled diagnostic '
        'tail fractions, not inferential p-values. A candidate entirely outside the initial 3-kb '
        'window is therefore not blocked by absence of a 3-kb signal.','',
        'Coding overlap, intervening genes, low-complexity repeats, a stronger <60-nt period, ambiguous '
        'array sequence, boundary saturation and initial/final copy-count disagreement prevent automatic '
        'acceptance and retain the candidate for review. Contig limitations elsewhere in the locus are '
        'reported without automatically rejecting an intact array. Two-copy candidates are explicitly '
        'exploratory, have no borrowed main-call p-value, and cannot be accepted.','',
        '## Alignment and full-array evaluation','',
        'Seed-centred segments are aligned with MAFFT before fitting the PWM. Conserved aligned bounds '
        'define a fixed-width PWM, calibrated conditionally against local-block shuffles for copy recovery. '
        'Copies are rescanned/chained and their common boundaries are checked again before exporting '
        'actual spacers. Gapped alignments are explicitly projected onto columns occupied by DNA in '
        'at least 80% of copies; minority-copy gaps are missing observations. The full alignment, '
        'projection mask and indel flags are retained. Such a fixed-width PWM may miss copies with '
        'internal indels, which can instead appear as missing copies or review flags. '
        'The complete-array information score is descriptive, not substituted '
        'into the pre-PWM seed-chain null. All coordinates are zero-based, half-open.','',
        '## Nulls and limitations','',
        'The local-block null preserves mononucleotide counts within successive 50-nt blocks. Exact '
        'dinucleotide shuffles preserve endpoints and every adjacent-base count. They use the Euler '
        'last-exit-tree construction, enumerating weighted arborescences over at most four DNA vertices. '
        'Tree multiplicity weights account for parallel edges. ACGT runs separated by ambiguity are '
        'shuffled separately with fixed endpoints. This is based on the Euler construction described by '
        '[Jiang et al., uShuffle](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/); it does not use '
        'an arbitrary random Euler walk. Both nulls rerun the entire seed-discovery/best-chain search.','',
        'Inference uses fixed-size Monte Carlo samples, avoiding an empirical p-value from the replication '
        'procedure\'s adaptive 200/2,000 stopping rule. Existing draws are reused only as a fixed prefix '
        'with verified input checksums; additional local draws follow the same seed schedule. The '
        'paper-style maximum-shuffle decision remains a distinct descriptive replication field. These '
        'are computational candidates, not biological validation or guarantees of precise RNA boundaries.','',
        'The exploratory file reports at most the configured number of highest-scoring two-copy tracts '
        'per locus; it is a prioritization list, not a complete calibrated two-copy census.','',
        '## Results','',f"Status counts: `{summary['assessment_counts']}`.",'',
        '| Locus | Copies | Adjusted p | Assessment | Flags |','|---|---:|---:|---|---|']
    lines.extend(f"| {r['locus_id']} | {r['copy_count']} | {float(r['bonferroni_p']):.4g} | {r['assessment_status']} | {r['flags'] or 'none'} |" for r in arrays)
    lines += ['', '## Files and reproduction','',
        '`locus_assessments.tsv`, `array_assessments.tsv`, `repeats.tsv`, `spacers.tsv`, the corresponding FASTAs, '
        '`alignments/`, `details/`, `background_tests.tsv`, `null_best_scores.tsv` and '
        '`exploratory_two_copy.tsv` contain the scores, coordinates, sequences, alignments and decisions. '
        '`accepted_arrays.fna`, `accepted_repeats.fna` and `accepted_spacers.fna` contain only automatically '
        'accepted candidates; other candidates remain in the complete assessment files.','',
        '```powershell','.venv/Scripts/python.exe scripts/art_array_assessment.py',
        '# Resume or redo downstream annotation using compatible saved nulls:',
        '.venv/Scripts/python.exe scripts/art_array_assessment.py --reuse-nulls','```','']
    (output/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,default=ROOT/'results/art_arrays')
    parser.add_argument('--loci',type=Path,default=ROOT/'results/oriented_art_loci')
    parser.add_argument('--output',type=Path,default=ROOT/'results/art_arrays_assessed')
    parser.add_argument('--config',type=Path,default=ROOT/'config/art_array_assessment.json')
    parser.add_argument('--reuse-nulls',action='store_true')
    run(parser.parse_args())


if __name__ == '__main__':
    main()
