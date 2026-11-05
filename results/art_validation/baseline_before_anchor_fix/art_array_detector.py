"""Custom ART DNA array scanner; internal coordinates are zero-based, half-open.

Implements independent published-style initial and extended discovery searches.
This is a documented reimplementation, not the authors' original code.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import numpy as np

from extract_art_loci import read_tsv, source_pieces, Segment, write_json, write_tsv

ROOT = Path(__file__).resolve().parents[1]
BASES = 'ACGT'
LUT = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(BASES):
    LUT[ord(b)] = i


def dna_array(sequence):
    return LUT[np.frombuffer(str(sequence).upper().encode('ascii'), dtype=np.uint8)]


def words(sequence, k):
    """Index valid ACGT words; ambiguous bases never match or seed a repeat."""
    seq = str(sequence).upper()
    index = defaultdict(list)
    code, valid = 0, 0
    mask = (1 << (2*k))-1
    for i, b in enumerate(seq):
        digit = BASES.find(b)
        if digit < 0:
            code, valid = 0, 0
        else:
            code = ((code << 2) | digit) & mask
            valid += 1
            if valid >= k:
                index[code].append(i-k+1)
    return index


def decode(code, k):
    return ''.join(BASES[(code >> (2*i)) & 3] for i in reversed(range(k)))


def nonoverlap(positions, k):
    result = []
    for p in sorted(set(positions)):
        if not result or p >= result[-1]+k:
            result.append(int(p))
    return result


def approximate_positions(index, code, k, mismatches):
    hits = list(index.get(code, ()))
    variants = []
    for i in range(k):
        old = (code >> (2*i)) & 3
        for digit in range(4):
            if digit != old:
                changed = code ^ ((old ^ digit) << (2*i))
                hits.extend(index.get(changed, ()))
                variants.append((i, changed))
    if mismatches == 2:
        for i, changed in variants:
            for j in range(i+1, k):
                old = (code >> (2*j)) & 3
                for digit in range(4):
                    if digit != old:
                        hits.extend(index.get(changed ^ ((old ^ digit) << (2*j)), ()))
    elif mismatches not in (0, 1):
        raise ValueError('Supported seed mismatch counts: 0, 1, 2')
    if mismatches == 0:
        hits = index.get(code, ())
    return nonoverlap(hits, k)


def best_chain(positions, low, high, tolerance=.30, missing=0, strict_gaps=False):
    """Longest chain under a common feasible spacing interval.

    Each edge constrains the same base spacing d by gap/(m*(1+t)) <= d
    <= gap/(m*(1-t)), m=1 or 2. Dynamic programming retains nondominated
    feasible intervals, so no arbitrary spacing grid is used. At most `missing`
    skipped copies are allowed over the entire chain.
    """
    pos = sorted(set(int(p) for p in positions))
    if len(pos) < 2:
        return dict(positions=pos, base_spacing=None, missing_copies=0)
    states = [[] for p in pos]
    best = ([], float(low), float(high), 0)
    for j, p in enumerate(pos):
        candidates = [(1, 0, float(low), float(high), (j,))]
        for i in range(j-1, -1, -1):
            gap = p-pos[i]
            if gap > high*(1+tolerance)*(2 if missing else 1)+1e-8:
                break
            if strict_gaps and not low <= gap <= high:
                continue
            for count, used, lo, hi, path in states[i]:
                for multiplier in (1, 2) if missing > used else (1,):
                    a = max(lo, gap/(multiplier*(1+tolerance)))
                    b = min(hi, gap/(multiplier*(1-tolerance)))
                    if a <= b+1e-8:
                        candidates.append((count+1, used+multiplier-1, a, b, path+(j,)))
        kept = []
        for state in sorted(candidates, key=lambda x: (-x[0], x[1], -(x[3]-x[2]), x[4])):
            count, used, lo, hi, path = state
            if any(c >= count and u <= used and a <= lo and b >= hi for c,u,a,b,_ in kept):
                continue
            kept.append(state)
            chain = [pos[i] for i in path]
            if (count, -used, -chain[0]) > (len(best[0]), -best[3], -(best[0][0] if best[0] else math.inf)):
                best = chain, lo, hi, used
        states[j] = kept
    chain, lo, hi, used = best
    # Pick the feasible spacing closest to the median observed (normalized) gap.
    gaps = np.diff(chain)
    middle = (lo+hi)/2
    normalized = gaps/np.where(gaps > middle*(1+tolerance), 2, 1)
    estimate = float(np.median(normalized)) if len(gaps) else middle
    spacing = min(hi, max(lo, estimate))
    return dict(positions=chain, base_spacing=spacing, missing_copies=used)


def dense_run(positions, minimum_gap, k):
    best = current = 1 if positions else 0
    for a, b in zip(positions, positions[1:]):
        current = current+1 if k <= b-a < minimum_gap else 1
        best = max(best, current)
    return best


def base_background(sequence):
    encoded = dna_array(sequence)
    counts = np.bincount(encoded[encoded < 4], minlength=4).astype(float) + .5
    return counts/counts.sum()


def conserved_block(matrix, anchor, fraction=.8, lapses=1):
    """Longest anchor-containing block with <=1 low-conservation column total."""
    if not len(matrix) or not matrix.shape[1]:
        return None
    counts = np.stack([(matrix == i).sum(axis=0) for i in range(4)], axis=1)
    good = counts.max(axis=1) >= fraction*len(matrix)-1e-9
    if not 0 <= anchor < len(good):
        return None
    prefix = np.r_[0, np.cumsum(~good)]
    best = None
    for a in np.flatnonzero(good[:anchor+1]):
        if prefix[anchor+1]-prefix[a] > lapses:
            continue
        b = min(len(good), int(np.searchsorted(prefix, prefix[a]+lapses, side='right'))-1)
        while b > anchor and not good[b-1]:
            b -= 1
        if b > anchor:
            key = (b-int(a), int(good[a:b].sum()), -int(a))
            if best is None or key > best[0]:
                best = key, int(a), int(b)
    if best is None:
        return None
    a, b = best[1:]
    return dict(start=a, end=b, counts=counts[a:b], conserved=good[a:b])


def boundary_estimate(sequence, positions, k, cfg):
    flank = cfg['boundary_flank_nt']
    encoded = dna_array(sequence)
    offsets = np.arange(-flank, k+flank)
    indices = np.array(positions)[:, None]+offsets
    valid = (indices >= 0) & (indices < len(encoded))
    matrix = np.full(indices.shape, 4, dtype=np.uint8)
    matrix[valid] = encoded[indices[valid]]
    block = conserved_block(matrix, flank+k//2, cfg['boundary_consensus_fraction'], cfg['boundary_allowed_lapses'])
    if block is None:
        return None
    a, b = block['start']-flank, block['end']-flank
    counts = block['counts']
    n = counts.sum(axis=1, keepdims=True)
    frequencies = np.divide(counts, n, out=np.zeros_like(counts, dtype=float), where=n > 0)
    bg = base_background(sequence)
    info = float(np.sum(np.where(frequencies > 0,
        frequencies*np.log2(np.maximum(frequencies, 1e-15)/bg), 0)))
    return dict(offset_start=a, offset_end=b, length_nt=b-a, information_bits=info,
                consensus=''.join(BASES[i] for i in counts.argmax(axis=1)),
                counts=counts.tolist(), boundary_at_search_limit=(block['start'] == 0 or block['end'] == len(offsets)),
                intervals=[[max(0, p+a), min(len(sequence), p+b)] for p in positions])


def initial_scan(sequence, cfg):
    k = cfg['initial_seed_nt']
    index = words(sequence, k)
    eligible = [(code, ps) for code, ps in index.items()
                if len(set(decode(code, k))) >= 3 and not re.search(r'(.)\1{5}', decode(code, k))]
    eligible.sort(key=lambda item: (-len(item[1]), item[0]))
    best = dict(R=0, positions=[], seed='', dense_run=0, dense_veto=False)
    max_dense = 0
    for code, ps in eligible[:cfg['initial_top_seeds']]:
        positions = approximate_positions(index, code, k, cfg['initial_max_mismatches'])
        chain = best_chain(positions, *cfg['initial_start_spacing_nt'],
                           cfg['initial_spacing_tolerance'], strict_gaps=True)
        dense = dense_run(positions, cfg['initial_start_spacing_nt'][0], k)
        max_dense = max(max_dense, dense)
        if (len(chain['positions']), len(ps), -code) > (best['R'], best.get('exact_seed_count', 0), -best.get('code', 0)):
            best = dict(R=len(chain['positions']), positions=chain['positions'], seed=decode(code, k),
                        dense_run=dense, code=code, exact_seed_count=len(ps), base_spacing=chain['base_spacing'])
    best['dense_veto'] = max_dense > best['R']
    best['maximum_dense_run_any_seed'] = max_dense
    best['seeds_tested'] = min(len(eligible), cfg['initial_top_seeds'])
    return best


def refinement_scan(sequence, cfg, return_candidates=True):
    k = cfg['refinement_seed_nt']
    index = words(sequence, k)
    candidates = []
    # "Recurring" explicitly means >=2 exact valid occurrences before
    # mismatch expansion. No low-complexity filter or top-N seed cap here.
    recurring = sorted(code for code, ps in index.items() if len(ps) >= 2)
    if len(set(str(sequence)) & set(BASES)) <= 1:
        return dict(score=0., candidates=[], seeds_tested=len(recurring))
    seen = set()
    for code in recurring:
        positions = approximate_positions(index, code, k, cfg['refinement_max_mismatches'])
        if len(positions) < cfg['minimum_copies_for_main_call']:
            continue
        chain = best_chain(positions, *cfg['refinement_base_spacing_nt'], cfg['refinement_spacing_tolerance'],
            missing=cfg['maximum_missing_copies_per_chain'] if cfg['allow_single_missing_copy'] else 0)
        selected = tuple(chain['positions'])
        if len(selected) < cfg['minimum_copies_for_main_call'] or selected in seen:
            continue
        seen.add(selected)
        boundary = boundary_estimate(sequence, selected, k, cfg)
        if boundary is None:
            continue
        score = (len(selected)-1)*boundary['information_bits']
        candidates.append(dict(score=score, seed=decode(code, k), **chain, boundary=boundary))
    candidates.sort(key=lambda c: (-c['score'], -len(c['positions']), c['positions'][0], c['seed']))
    return dict(score=candidates[0]['score'] if candidates else 0.,
                candidates=candidates if return_candidates else [], seeds_tested=len(recurring))


def shuffled(sequence, rng, block=None):
    data = np.array(list(sequence), dtype='U1')
    if block is None:
        rng.shuffle(data)
    else:
        for start in range(0, len(data), block):
            rng.shuffle(data[start:start+block])
    return ''.join(data)


def null_batch(task):
    sequence, mode, cfg, seeds = task
    result = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        seq = shuffled(sequence, rng, None if mode == 'initial' else cfg['shuffle_block_nt'])
        score = initial_scan(seq, cfg)['R'] if mode == 'initial' else refinement_scan(seq, cfg, False)['score']
        result.append((seed, score))
    return result


def stable_seed(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(map(str, parts)).encode()).digest()[:8], 'little')


def null_scores(sequence, mode, cfg, n, locus, window, pool, start=0):
    seeds = [stable_seed(cfg['random_seed'], locus, mode, window, i) for i in range(start, n)]
    tasks = [(sequence, mode, cfg, seeds[i:i+20]) for i in range(0, len(seeds), 20)]
    rows = []
    for batch in pool.map(null_batch, tasks):
        rows.extend(batch)
    return rows


def empirical(score, scores):
    exceedances = sum(s >= score for s in scores)
    return dict(observed_score=score, null_max=max(scores, default=0), shuffles=len(scores),
                exceedances=exceedances, empirical_p=(exceedances+1)/(len(scores)+1),
                exceeds_all_shuffles=bool(scores) and score > max(scores))


def context_flags(intervals, features, rt_start, cfg):
    first, last = intervals[0][0], intervals[-1][1]
    starts = [a for a,b in intervals]
    coding = []
    for f in features:
        if f.type == 'CDS' and all(any(int(p.start) <= a and b <= int(p.end) for p in f.location.parts)
                                   for a,b in intervals) and all(int(g) % 3 == 0 for g in np.diff(starts)):
            coding.append(f.qualifiers.get('protein_id', f.qualifiers.get('locus_tag', ['unknown']))[0])
    intervening = []
    for f in features:
        if f.type == 'CDS' and len(f.location) >= cfg['adjacency_gene_min_nt'] and (
                int(f.location.start) >= last and int(f.location.end) <= rt_start):
            intervening.append(f.qualifiers.get('protein_id', f.qualifiers.get('locus_tag', ['unknown']))[0])
    return dict(coding_repeat_gene_ids=sorted(set(coding)), intervening_gene_ids=sorted(set(intervening)),
                coding_repeat=bool(coding), adjacent_to_RT=not intervening,
                array_start0=first, array_end0=last, array_span_nt=last-first,
                distance_to_RT_nt=rt_start-last)


def pwm_scores(sequence, pwm):
    encoded = dna_array(sequence)
    width = len(pwm)
    n = len(encoded)-width+1
    if n <= 0:
        return np.array([])
    result = np.zeros(n)
    invalid = np.zeros(n, dtype=bool)
    for j in range(width):
        col = encoded[j:j+n]
        invalid |= col >= 4
        result += pwm[j, np.minimum(col, 3)]
    result[invalid] = -np.inf
    return result


def nonoverlap_scored(positions, scores, width):
    selected = []
    for p in sorted(positions, key=lambda p: (-scores[p], p)):
        if all(p+width <= old or old+width <= p for old in selected):
            selected.append(int(p))
    return sorted(selected)


def pwm_refine(sequence, candidate, cfg, locus):
    counts = np.array(candidate['boundary']['counts'], dtype=float)
    counts += cfg['pwm_pseudocount_per_base']
    probabilities = counts/counts.sum(axis=1, keepdims=True)
    pwm = np.log2(probabilities/base_background(sequence))
    null = []
    for i in range(cfg['pwm_shuffles']):
        rng = np.random.default_rng(stable_seed(cfg['random_seed'], locus, 'PWM', candidate['seed'], i))
        values = pwm_scores(shuffled(sequence, rng, cfg['shuffle_block_nt']), pwm)
        null.append(float(values.max()) if len(values) else -math.inf)
    scores = pwm_scores(sequence, pwm)
    threshold = max(null)
    hits = nonoverlap_scored(np.flatnonzero(scores > threshold), scores, len(pwm))
    chain = best_chain(hits, *cfg['refinement_base_spacing_nt'], cfg['refinement_spacing_tolerance'],
        missing=cfg['maximum_missing_copies_per_chain'] if cfg['allow_single_missing_copy'] else 0)
    return dict(pwm=pwm.tolist(), null_max=threshold, null_max_scores=null,
                significant_hit_starts=hits, periodic_hit_starts=chain['positions'],
                off_chain_significant_starts=sorted(set(hits)-set(chain['positions'])),
                hit_scores=[float(scores[p]) for p in hits], chain=chain)


def align_boundaries(sequence, seed_positions, k, cfg, output, array_id):
    flank = cfg['boundary_flank_nt']
    fragments, starts = [], []
    for i, p in enumerate(seed_positions):
        a, b = max(0, p-flank), min(len(sequence), p+k+flank)
        starts.append(a)
        fragments.append(SeqRecord(Seq(sequence[a:b]), id=f'copy{i+1}', description=''))
    raw = output/'alignments'/f'{array_id}.flanks.fna'
    aligned = output/'alignments'/f'{array_id}.mafft.fna'
    SeqIO.write(fragments, raw, 'fasta')
    command = ['cmd.exe', '/d', '/c', str(ROOT/'.tools/mafft-win/mafft.bat'),
               '--auto', '--thread', '1', '--inputorder', '--out', str(aligned), str(raw)]
    with (output/'alignments'/f'{array_id}.mafft.log').open('w', encoding='utf-8') as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    with aligned.open(encoding='utf-8') as handle:
        records = list(SeqIO.parse(handle, 'fasta'))
    by_id = {r.id: str(r.seq).upper() for r in records}
    texts = [by_id[r.id] for r in fragments]
    matrix = np.array([[BASES.find(c) if c in BASES else 4 for c in text] for text in texts], dtype=np.uint8)
    coordinates = []
    for text, start in zip(texts, starts):
        row, cursor = [], start
        for char in text:
            row.append(cursor if char != '-' else None)
            cursor += char != '-'
        coordinates.append(row)
    anchor = next(i for i,p in enumerate(coordinates[0]) if p == seed_positions[0]+k//2)
    block = conserved_block(matrix, anchor, cfg['boundary_consensus_fraction'], cfg['boundary_allowed_lapses'])
    if block is None:
        return dict(status='no_conserved_block_after_MAFFT', command=command)
    a, b = block['start'], block['end']
    intervals = []
    aligned_repeats = []
    for text, coords in zip(texts, coordinates):
        positions = [p for p in coords[a:b] if p is not None]
        if not positions:
            return dict(status='empty_copy_after_MAFFT', command=command)
        intervals.append([min(positions), max(positions)+1])
        aligned_repeats.append(text[a:b])
    return dict(status='aligned', command=command, intervals=intervals,
        aligned_repeat_sequences=aligned_repeats, alignment_columns=b-a,
        lengths_nt=[y-x for x,y in intervals],
        consensus=''.join(BASES[i] for i in block['counts'].argmax(axis=1)),
        boundary_at_search_limit=(a == 0 or b == matrix.shape[1]),
        gaps_in_repeat=any('-' in text for text in aligned_repeats))


def distinct_arrays(candidates):
    """Merge seeds describing the same tract; prefer more copies, then score."""
    result = []
    for candidate in sorted(candidates, key=lambda c: (-len(c['positions']), -c['score'], c['positions'])):
        intervals = candidate['boundary']['intervals']
        if any(sum(any(max(a,x) < min(b,y) for x,y in old['boundary']['intervals'])
                   for a,b in intervals) >= 2 for old in result):
            continue
        result.append(candidate)
    return result


def initial_support(initial, passed, initial_offset, intervals, k):
    if not passed:
        return False
    positions = [p+initial_offset for p in initial['positions']]
    return sum(any(max(p,a) < min(p+k,b) for a,b in intervals) for p in positions) >= 2


def run(args):
    started = time.perf_counter()
    cfg = json.loads(args.config.read_text(encoding='utf-8'))['art_array_search']
    if cfg['enforce_observed_ranges_as_filters']:
        raise ValueError('Observed ranges are descriptive in this implementation; hard filters are disabled')
    if not cfg['estimate_repeat_boundaries']:
        raise ValueError('Full repeat boundary estimation is required to measure spacers')
    output = args.output.resolve()
    for folder in (output, output/'alignments', output/'loci', output/'code', output/'details'):
        folder.mkdir(parents=True, exist_ok=True)
    inputs = args.loci.resolve()
    manifests = read_tsv(inputs/'locus_manifest.tsv')
    raw_maps = read_tsv(inputs/'coordinate_map.tsv')
    maps = defaultdict(list)
    for row in raw_maps:
        maps[row['locus_id']].append(Segment(**{k: int(row[k]) for k in (
            'source_start0','source_end0','locus_start0','locus_end0','orientation')}))
    for segments in maps.values():
        segments.sort(key=lambda s: s.locus_start0)
    locus_rows, array_rows, repeat_rows, spacer_rows, null_rows, candidate_rows = [], [], [], [], [], []
    repeat_fasta, spacer_fasta = [], []
    commands = []
    with ProcessPoolExecutor(max_workers=cfg['workers']) as pool:
        for number, meta in enumerate(manifests, 1):
            locus_id = meta['locus_id']
            record = SeqIO.read(inputs/meta['genbank_path'], 'genbank')
            rt_start = int(meta['RT_locus_start0'])
            ext_start = max(0, rt_start-cfg['extended_window_nt'])
            ini_start = max(0, rt_start-cfg['initial_window_nt'])
            initial_sequence = str(record.seq[ini_start:rt_start]).upper()
            extended_sequence = str(record.seq[ext_start:rt_start]).upper()
            print(f'[{number}/{len(manifests)}] {locus_id}: independent initial and extended scans', flush=True)
            initial = initial_scan(initial_sequence, cfg)
            refinement = refinement_scan(extended_sequence, cfg)
            ref3 = refinement_scan(initial_sequence, cfg)
            initial_null = null_scores(initial_sequence, 'initial', cfg, cfg['initial_shuffles'],
                                       locus_id, len(initial_sequence), pool)
            initial_test = empirical(initial['R'], [s for _,s in initial_null])
            initial_pass = (initial['R'] >= cfg['minimum_copies_for_main_call'] and not initial['dense_veto']
                            and initial_test['exceeds_all_shuffles'])
            for i,(seed,score) in enumerate(initial_null):
                null_rows.append(dict(locus_id=locus_id, stage='initial_mononucleotide', window_nt=len(initial_sequence),
                                      shuffle_index=i, rng_seed=seed, score=score))
            nulls = {}
            for label, sequence in (('3000', initial_sequence), ('6000', extended_sequence)):
                nulls[label] = null_scores(sequence, 'refinement_'+label, cfg, cfg['initial_refinement_shuffles'],
                                           locus_id, len(sequence), pool)
            preliminary_max = max((score for rows in nulls.values() for seed,score in rows), default=0)
            # A failed strict maximum comparison cannot be rescued by adding more
            # shuffles; only provisionally passing borderline chains need retest.
            borderline = any(preliminary_max < c['score'] < cfg['borderline_score_ratio']*preliminary_max
                             for c in refinement['candidates'])
            if borderline:
                print(f'  borderline chain: increasing each refinement null to {cfg["borderline_refinement_shuffles"]}', flush=True)
                for label, sequence in (('3000', initial_sequence), ('6000', extended_sequence)):
                    nulls[label] += null_scores(sequence, 'refinement_'+label, cfg, cfg['borderline_refinement_shuffles'],
                        locus_id, len(sequence), pool, start=cfg['initial_refinement_shuffles'])
            for label, rows in nulls.items():
                for i,(seed,score) in enumerate(rows):
                    null_rows.append(dict(locus_id=locus_id, stage='refinement_local_blocks_'+label,
                        window_nt=len(initial_sequence if label == '3000' else extended_sequence),
                        shuffle_index=i, rng_seed=seed, score=score))
            null_max = max((score for rows in nulls.values() for seed,score in rows), default=0)
            candidates = []
            for candidate in refinement['candidates']:
                intervals = [[a+ext_start,b+ext_start] for a,b in candidate['boundary']['intervals']]
                flags = context_flags(intervals, record.features, rt_start, cfg)
                passed = candidate['score'] > null_max
                candidate_rows.append(dict(locus_id=locus_id, seed=candidate['seed'], copies=len(candidate['positions']),
                    score=candidate['score'], exceeds_both_null_maxima=passed, **flags,
                    seed_positions0_json=json.dumps([p+ext_start for p in candidate['positions']]),
                    boundary_length_nt=candidate['boundary']['length_nt']))
                if passed and not flags['coding_repeat'] and flags['adjacent_to_RT']:
                    candidates.append(candidate)
            accepted = []
            pwm_details = []
            # Deduplicate before PWM construction, so shifted seeds do not each
            # fit an extra repeat model. All candidate scores remain available.
            for candidate_index, candidate in enumerate(distinct_arrays(candidates), 1):
                # Failed models keep their own identity in null/provenance rows.
                array_id = f'{locus_id}__array{candidate_index}'
                pwm = pwm_refine(extended_sequence, candidate, cfg, locus_id)
                pwm_details.append(dict(array_id=array_id, seed=candidate['seed'], **pwm))
                for i,score in enumerate(pwm['null_max_scores']):
                    null_rows.append(dict(locus_id=locus_id, stage='PWM_'+array_id,
                        window_nt=len(extended_sequence), shuffle_index=i,
                        rng_seed=stable_seed(cfg['random_seed'], locus_id, 'PWM', candidate['seed'], i), score=score))
                if len(pwm['periodic_hit_starts']) < cfg['minimum_copies_for_main_call']:
                    continue
                seed_positions = [p-candidate['boundary']['offset_start'] for p in pwm['periodic_hit_starts']]
                ungapped = boundary_estimate(extended_sequence, seed_positions, cfg['refinement_seed_nt'], cfg)
                if ungapped is None:
                    continue
                aligned = align_boundaries(extended_sequence, seed_positions, cfg['refinement_seed_nt'], cfg, output, array_id)
                commands.append(aligned['command'])
                if aligned['status'] != 'aligned':
                    continue
                intervals = [[a+ext_start, b+ext_start] for a,b in aligned['intervals']]
                if any(b > c for (a,b),(c,d) in zip(intervals, intervals[1:])):
                    continue
                flags = context_flags(intervals, record.features, rt_start, cfg)
                if flags['coding_repeat'] or not flags['adjacent_to_RT']:
                    continue
                spacings = [c-a for (a,b),(c,d) in zip(intervals, intervals[1:])]
                spacers = [c-b for (a,b),(c,d) in zip(intervals, intervals[1:])]
                repeat_lengths = [b-a for a,b in intervals]
                supported_initial = initial_support(initial, initial_pass, ini_start, intervals, cfg['initial_seed_nt'])
                mode_label = 'initial_mode_and_extended' if supported_initial else 'extended_only'
                tests = {label: empirical(candidate['score'], [s for _,s in rows]) for label,rows in nulls.items()}
                row = dict(array_id=array_id, locus_id=locus_id, accession=meta['accession'], RT_protein_id=meta['RT_protein_id'],
                    discovery_label=mode_label, initial_mode_pass=supported_initial, locus_initial_mode_pass=initial_pass,
                    copy_count=len(intervals), **flags, chain_score=candidate['score'],
                    seed=candidate['seed'], seed_chain_copies=len(candidate['positions']),
                    seed_missing_copies=candidate['missing_copies'],
                    repeat_consensus=aligned['consensus'], repeat_alignment_columns=aligned['alignment_columns'],
                    repeat_lengths_nt_json=json.dumps(repeat_lengths),
                    ungapped_repeat_length_nt=ungapped['length_nt'],
                    repeat_start_spacings_nt_json=json.dumps(spacings), spacer_lengths_nt_json=json.dumps(spacers),
                    repeat_length_outside_observed=any(not cfg['observed_repeat_length_nt'][0] <= x <= cfg['observed_repeat_length_nt'][1]
                                                     for x in repeat_lengths),
                    spacer_length_outside_observed=any(not cfg['observed_spacer_length_nt'][0] <= x <= cfg['observed_spacer_length_nt'][1]
                                                     for x in spacers),
                    boundary_at_search_limit=aligned['boundary_at_search_limit'] or ungapped['boundary_at_search_limit'],
                    gapped_repeat_alignment=aligned['gaps_in_repeat'],
                    null_3000_max=tests['3000']['null_max'], null_6000_max=tests['6000']['null_max'],
                    null_shuffles_per_window=len(nulls['6000']),
                    empirical_p_3000=tests['3000']['empirical_p'], empirical_p_6000=tests['6000']['empirical_p'],
                    PWM_null_max=pwm['null_max'], PWM_significant_off_chain=len(pwm['off_chain_significant_starts']))
                accepted.append(row)
                array_rows.append(row)
                write_json(output/'details'/f'{array_id}.json', dict(summary=row, seed_chain=candidate,
                    final_ungapped=ungapped, final_MAFFT=aligned, pwm=pwm, shuffle_tests=tests))
                for i,((a,b),aligned_text) in enumerate(zip(intervals, aligned['aligned_repeat_sequences']), 1):
                    repeat_id = f'{array_id}.repeat{i}'
                    text = str(record.seq[a:b])
                    assert text == aligned_text.replace('-', '')
                    original = source_pieces(a,b,maps[locus_id])
                    assert sum(p['source_end0']-p['source_start0'] for p in original) == b-a
                    repeat_rows.append(dict(repeat_id=repeat_id, array_id=array_id, locus_id=locus_id, copy_index=i,
                        locus_start0=a, locus_end0=b, length_nt=b-a, sequence=text,
                        aligned_sequence=aligned_text, source_segments_json=json.dumps(original)))
                    repeat_fasta.append(SeqRecord(Seq(text), id=repeat_id, description=f'locus=[{a},{b})'))
                    record.features.append(SeqFeature(SimpleLocation(a,b,strand=1), type='repeat_region',
                        qualifiers=dict(note=['Custom ART array candidate repeat; not a functional RNA annotation'],
                                        rpt_family=[array_id], rpt_type=['direct'])))
                    if i < len(intervals):
                        c = intervals[i][0]
                        spacer_id = f'{array_id}.spacer{i}'
                        spacer_text = str(record.seq[b:c])
                        spacer_rows.append(dict(spacer_id=spacer_id, array_id=array_id, locus_id=locus_id,
                            spacer_index=i, locus_start0=b, locus_end0=c, length_nt=c-b,
                            repeat_start_spacing_nt=c-a, preceding_repeat_length_nt=b-a,
                            source_segments_json=json.dumps(source_pieces(b,c,maps[locus_id])), sequence=spacer_text))
                        assert c-b == (c-a)-(b-a)
                        spacer_fasta.append(SeqRecord(Seq(spacer_text), id=spacer_id, description=f'locus=[{b},{c})'))
            SeqIO.write(record, output/'loci'/f'{locus_id}.arrays.gb', 'genbank')
            status = ('array_candidate_detected' if accepted else
                      'initial_only_not_confirmed_by_delimitation' if initial_pass else
                      'not_assessed_short_upstream' if rt_start < cfg['minimum_assessable_upstream_nt'] else
                      'no_supported_array_detected_not_evidence_of_absence')
            locus_rows.append(dict(locus_id=locus_id, RT_protein_id=meta['RT_protein_id'],
                upstream_available_nt=rt_start, initial_window_nt=len(initial_sequence), extended_window_nt=len(extended_sequence),
                extended_search_executed=True, initial_R=initial['R'], initial_seed=initial['seed'],
                initial_dense_veto=initial['dense_veto'], initial_null_max=initial_test['null_max'],
                initial_empirical_p=initial_test['empirical_p'], initial_pass=initial_pass,
                refinement_3000_best_score=ref3['score'], refinement_6000_best_score=refinement['score'],
                refinement_3000_seeds=ref3['seeds_tested'], refinement_6000_seeds=refinement['seeds_tested'],
                refinement_null_3000_max=max(s for _,s in nulls['3000']),
                refinement_null_6000_max=max(s for _,s in nulls['6000']),
                refinement_shuffles_per_window=len(nulls['6000']), borderline_retested=borderline,
                array_count=len(accepted), extended_only_array_count=sum(r['discovery_label'] == 'extended_only' for r in accepted),
                status=status))
            write_json(output/'details'/f'{locus_id}.scan.json', dict(initial=initial, initial_test=initial_test,
                refinement_3000=ref3, refinement_6000=refinement, PWM_models=pwm_details,
                extended_search_independent_of_initial=True))
            print(f'  initial pass={initial_pass}, delimited arrays={len(accepted)}, status={status}', flush=True)
    for name, rows, empty in (
        ('locus_calls.tsv', locus_rows, ['locus_id','status']),
        ('arrays.tsv', array_rows, ['array_id','locus_id','copy_count']),
        ('repeats.tsv', repeat_rows, ['repeat_id','array_id','locus_start0','locus_end0']),
        ('spacers.tsv', spacer_rows, ['spacer_id','array_id','locus_start0','locus_end0']),
        ('null_scores.tsv', null_rows, ['locus_id','stage','score']),
        ('all_seed_chains.tsv', candidate_rows, ['locus_id','seed','score'])):
        write_tsv(output/name, rows, list(rows[0]) if rows else empty)
    SeqIO.write(repeat_fasta, output/'repeats.fna', 'fasta')
    SeqIO.write(spacer_fasta, output/'spacers.fna', 'fasta')
    write_json(output/'run_config.json', {'art_array_search': cfg})
    write_json(output/'commands.json', commands)
    for script in (Path(__file__), Path(__file__).with_name('test_art_array_detector.py')):
        if script.exists():
            shutil.copyfile(script, output/'code'/script.name)
    # The coordinate helper is part of the complete reproducible source snapshot.
    shutil.copyfile(Path(__file__).with_name('extract_art_loci.py'), output/'code/extract_art_loci.py')
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(), python=sys.version,
        biopython=Bio.__version__, numpy=np.__version__, random_generator='numpy.PCG64', command=sys.argv,
        paper_path='literature/Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf',
        paper_sha256=hashlib.sha256((ROOT/'literature/Yoon et al. 2026 Array Associated Reverse Transcriptases (2).pdf').read_bytes()).hexdigest(),
        input_hashes={p.relative_to(inputs).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [inputs/'locus_manifest.tsv', inputs/'coordinate_map.tsv']+
                      [inputs/m['genbank_path'] for m in manifests]},
        mafft_launcher_sha256=hashlib.sha256((ROOT/'.tools/mafft-win/mafft.bat').read_bytes()).hexdigest())
    write_json(output/'provenance.json', provenance)
    summary = dict(loci_scanned=len(locus_rows), initial_mode_pass=sum(r['initial_pass'] for r in locus_rows),
        delimited_arrays=len(array_rows), loci_with_delimited_arrays=sum(r['array_count'] > 0 for r in locus_rows),
        extended_only_arrays=sum(r['discovery_label'] == 'extended_only' for r in array_rows),
        repeat_copies=len(repeat_rows), spacers=len(spacer_rows),
        borderline_loci_retested=sum(r['borderline_retested'] for r in locus_rows),
        shuffled_scans=len(null_rows), elapsed_seconds=round(time.perf_counter()-started, 2),
        coordinate_convention='zero_based_half_open',
        interpretation='computational_ART_array_candidates_not_experimental_or_functional_validation')
    write_json(output/'summary.json', summary)
    write_json(output/'verification.json', dict(passed=True, all_loci_extended_searched=True,
        all_called_arrays_exceed_both_refinement_null_maxima=all(r['chain_score'] > max(r['null_3000_max'],r['null_6000_max']) for r in array_rows),
        all_main_calls_have_at_least_three_copies=all(r['copy_count'] >= cfg['minimum_copies_for_main_call'] for r in array_rows),
        repeat_and_spacer_sequence_coordinates_checked=True, source_mapping_lengths_checked=True))
    write_report(output, summary, locus_rows, array_rows, cfg)
    from art_array_initial_candidates import export_initial_candidates
    summary = export_initial_candidates(output, inputs, cfg)
    summary['elapsed_seconds'] = round(time.perf_counter()-started, 2)
    write_json(output/'summary.json', summary)
    write_json(output/'output_sha256.json', {p.relative_to(output).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'output_sha256.json'})
    print(json.dumps(summary, indent=2), flush=True)


def write_report(output, summary, loci, arrays, cfg):
    lines = ['# Dedicated ART array detector', '',
        'Custom DNA repeat scanning based on Yoon et al., Methods pp. 32-33. This is a documented '
        'reimplementation; the authors\' executable code and underspecified choices were not available.', '',
        'Both modes run independently on every RT-oriented upstream locus. A failed 3-kb initial scan '
        'never prevents the 6-kb scan. `initial_mode_and_extended` means both modes passed this '
        'implementation; `extended_only` means only the extended procedure yielded a delimited candidate. '
        'Neither label claims exact numerical reproduction of the authors\' analysis.', '',
        '## Results', '',
        f"Scanned {summary['loci_scanned']} loci; initial-mode passes: {summary['initial_mode_pass']}; "
        f"delimited arrays: {summary['delimited_arrays']}; extended-only arrays: {summary['extended_only_arrays']}. "
        f"Retested {summary['borderline_loci_retested']} borderline loci with "
        f"{cfg['borderline_refinement_shuffles']} shuffles per window.", '',
        '| Locus | Initial pass | Delimited arrays | Status |', '|---|---|---:|---|']
    lines += [f"| {r['locus_id']} | {r['initial_pass']} | {r['array_count']} | {r['status']} |" for r in loci]
    lines += ['', '## Outputs', '',
        '- `locus_calls.tsv`: both modes, failures, short upstream flanks and null maxima.',
        '- `arrays.tsv`: final boundaries, copy counts, repeat lengths, start spacings, actual spacer lengths and range flags.',
        '- `repeats.tsv`, `spacers.tsv`, and their FASTAs: sequences and original-source coordinate segments.',
        '- `all_seed_chains.tsv`: every scored extended candidate, including coding/non-adjacent/rejected chains.',
        '- `null_scores.tsv`: every initial, refinement and PWM null maximum, with reproducible RNG seeds.',
        '- `alignments/`: ungapped flanks, MAFFT alignments and logs; `details/`: seed, PWM and boundary evidence.',
        '- `loci/`: original oriented annotations plus candidate repeat_region features.', '',
        'All TSV/JSON positions are zero-based, half-open in the oriented locus; source segments explicitly '
        'record the orientation back to the original contig. GenBank files use one-based inclusive locations. '
        'Spacer length = next repeat start - current full repeat end. Start spacing is reported separately.', '',
        '## Implementation choices and limits', '',
        '- Initial seeds: top 20 exact-frequency 14-mers after the published diversity/homopolymer filters; '
        'up to 2 mismatches, non-overlapping matches, gaps restricted to 100-450 nt. The common spacing '
        'tolerance is 30% (the initial-method paragraph does not specify a numerical tolerance). '
        'The dense veto compares the longest <100-nt run across all tested seeds with the best regular run. '
        'The exceptional one-off exact-12-mer rescue in the paper is not applied automatically.',
        '- Extended seeds: every 10-mer with >=2 exact valid occurrences before mismatch expansion, '
        'up to 1 mismatch, no top-seed cap or initial-mode low-complexity filter. All recurring seeds '
        'are rescanned in every shuffled window. Ambiguous bases never match. Exact recurrence as '
        'the seed eligibility rule is explicit because the paper does not specify that detail.',
        '- Chains use a shared feasible spacing interval, a base spacing of 60-600 nt with 30% tolerance, '
        'and at most one skipped copy over a whole chain. This is a conservative interpretation of '
        '"single skipped copies allowed". Longest chains are retained per seed; their score is '
        '(copies-1) times summed column KL information relative to regional mononucleotide background. '
        'The best score across all seeds defines each shuffle maximum. The paper does not give '
        'the information-content formula or tie-breaking rules.',
        '- Conserved boundaries are the longest block containing the seed midpoint with >=80% consensus '
        'per column, allowing one low-conservation column, with conserved endpoints. Preliminary '
        f"ungapped flanks extend {cfg['boundary_flank_nt']} nt either side of each seed. A PWM "
        'with 0.5 pseudocount per base is calibrated against 200 new local-block shuffles; above-maximum '
        'non-overlapping matches are chained again. The final boundaries are measured in ungapped '
        'copies and MAFFT-aligned flanks; aligned bounds determine exported repeats/spacers. Boundary '
        'search-limit flags are recorded. This finite flank is a computational search parameter, not '
        'an asserted biological maximum repeat length.',
        '- The initial statistic exceeds all 100 mononucleotide shuffle run lengths. Extended candidates '
        'must strictly exceed local-block shuffle maxima in BOTH the 3-kb and 6-kb backgrounds, initially '
        '200 each. Provisionally passing chains below twice the initial maximum trigger extension '
        'to 2,000 total shuffles per window (the first 200 are retained). Already failed chains cannot '
        'be rescued by adding shuffles and are not retested. Independent 50-nt blocks retain their '
        'own base composition, including the final short block.',
        '- Empirical p = (1 + null scores >= observed)/(1 + shuffle count). These are finite-null '
        'comparisons for each window, not probabilities of ART identity or dataset-wide FDR. Null '
        'search maxima account for the within-window seed search. PWM thresholds are conditional on '
        'the fitted repeat model and are not an independent biological validation.',
        '- Chains inside one annotated CDS with all start spacings divisible by three are set aside '
        'as coding repeats. A CDS >=300 nt wholly between the final repeat and RT prevents adjacency. '
        'Shifted seed descriptions of a tract are merged when at least two repeat intervals overlap; '
        'the chain with more copies, then higher score, supplies the PWM. Full seed evidence is retained.',
        '- Published 15-49 nt repeats and 120-220 nt spacers are descriptive flags, never hard filters. '
        'Array spans are not capped at 4.1 kb. Failure to detect is not evidence of biological absence. '
        'No folding, RNA processing, protein binding or functional activity is inferred.', '',
        '## Reproduction', '',
        '```powershell', '.venv/Scripts/python.exe scripts/art_array_detector.py',
        '.venv/Scripts/python.exe -m unittest discover -s scripts -p test_art_array_detector.py -v', '```', '']
    (output/'REPORT.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT/'config/art_array_search.json')
    parser.add_argument('--loci', type=Path, default=ROOT/'results/oriented_art_loci')
    parser.add_argument('--output', type=Path, default=ROOT/'results/art_arrays')
    parser.add_argument('--finalize-only', action='store_true',
                        help='Export initial boundaries/review flags from an unchanged, completed calibrated run')
    args = parser.parse_args()
    if args.finalize_only:
        output, inputs = args.output.resolve(), args.loci.resolve()
        cfg = json.loads(args.config.read_text())['art_array_search']
        if json.loads((output/'run_config.json').read_text())['art_array_search'] != cfg:
            raise ValueError('Finalization cannot change calibrated configuration')
        provenance = json.loads((output/'provenance.json').read_text())
        for path, digest in provenance['input_hashes'].items():
            if hashlib.sha256((inputs/path).read_bytes()).hexdigest() != digest:
                raise ValueError(f'Calibrated input changed: {path}')
        previous_hashes = json.loads((output/'output_sha256.json').read_text())
        for relative in ('arrays.tsv','repeats.tsv','spacers.tsv','null_scores.tsv','locus_calls.tsv'):
            if hashlib.sha256((output/relative).read_bytes()).hexdigest() != previous_hashes[relative]:
                raise ValueError(f'Calibrated output changed: {relative}')
        summary = json.loads((output/'summary.json').read_text())
        write_report(output,summary,read_tsv(output/'locus_calls.tsv'),read_tsv(output/'arrays.tsv'),cfg)
        from art_array_initial_candidates import export_initial_candidates
        summary = export_initial_candidates(output,inputs,cfg)
        provenance['finalization_command'] = sys.argv
        provenance['finalized_utc'] = datetime.now(timezone.utc).isoformat()
        write_json(output/'provenance.json',provenance)
        for script in (Path(__file__),Path(__file__).with_name('test_art_array_detector.py')):
            shutil.copyfile(script,output/'code'/script.name)
        write_json(output/'output_sha256.json', {p.relative_to(output).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'output_sha256.json'})
        print(json.dumps(summary,indent=2))
    else:
        run(args)


if __name__ == '__main__':
    main()
