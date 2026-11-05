"""Tuned PILER-CR inside the same selection-aware shuffle wrapper as detector D.

For the 31 baseline-comparison inputs with copy-level truth (24 designed cases with
planted repeat intervals and the seven published reference loci), regenerates
exactly the seeded local-block and dinucleotide shuffles whose D scores are stored,
verifies the regeneration by rescoring a subset with D, runs tuned PILER-CR on the
observed window and every shuffle, and computes the same conservative empirical
p = max(p_local, p_dinucleotide) with PILER-CR's best-array score as the statistic.
PILER-CR 1.06 segfaults deterministically on a few inputs; a crashed shuffle counts as
an exceedance (conservative, primary) or as no array (sensitivity), and both are reported.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from Bio import SeqIO

from art_array_detector import ROOT, refinement_scan, shuffled, empirical
from art_array_assessment import dinucleotide_shuffle
from array_finder_tools import run_pilercr, best_array_score, PILERCR_TUNED
from extract_art_loci import read_tsv, write_tsv, write_json

OUT = ROOT/'results/array_calibration_20261105/pilercr_shuffle_wrapper'
VALIDATION = ROOT/'results/art_validation'
PLANTED = ('long_spacer_120', 'long_spacer_180', 'long_spacer_220', 'imperfect_copies', 'missing_copy',
           'extended_beyond_3kb', 'unrelated_RT_with_array', 'repetitive_DNA')
KINDS = ('local_blocks', 'dinucleotide')
VERIFY_DRAWS = (0, 1, 199, 200, 500, 999)


def search_config():
    return json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']


def inputs():
    """(input_id, group, 6-kb window, {kind: [(seed, stored D score)]}, stored D observed score)."""
    rows = []
    cases = {r['case_id']: r for r in read_tsv(VALIDATION/'sequence_cases.tsv')}
    for case, row in cases.items():
        if row['test'] not in PLANTED:
            continue
        record = SeqIO.read(VALIDATION/'loci'/(case+'.gb'), 'genbank')
        rt = next(int(f.location.start) for f in record.features if f.type == 'CDS' and f.qualifiers.get('protein_id') == ['synthetic_RT'])
        sequence = str(record.seq[max(0, rt-6000):rt]).upper()
        cache = json.loads((VALIDATION/'nulls'/(case+'.json')).read_text())
        if cache['key']['sequence_sha256'] != hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError('Window does not match its stored nulls: '+case)
        rows.append((case, 'designed_'+row['test'], sequence, {k: cache['scores'][k] for k in KINDS}, float(row['best_search_score'])))
    manifest = {r['locus_id']: r for r in read_tsv(ROOT/'results/oriented_art_loci/locus_manifest.tsv')}
    tests = {r['locus_id']: r for r in read_tsv(ROOT/'results/art_arrays_assessed/background_tests.tsv')}
    for ref in read_tsv(VALIDATION/'reference_recovery.tsv'):
        if ref['gold_status'] != 'published_core_anchored':
            continue
        locus = ref['locus_id']
        record = SeqIO.read(ROOT/'results/oriented_art_loci'/manifest[locus]['genbank_path'], 'genbank')
        rt = int(manifest[locus]['RT_locus_start0'])
        sequence = str(record.seq[max(0, rt-6000):rt]).upper()
        cache = json.loads((ROOT/'results/art_arrays_assessed/nulls'/(locus+'.json')).read_text())
        if cache['key']['sequence_sha256']['6000'] != hashlib.sha256(sequence.encode()).hexdigest():
            raise ValueError('Window does not match its stored nulls: '+locus)
        rows.append((locus, 'published_reference', sequence, {k: cache['scores'][k+'_6000'] for k in KINDS},
                     float(tests[locus]['local_blocks_6000_best_score'])))
    if len(rows) != 31:
        raise ValueError(f'Expected 31 inputs with copy-level truth, found {len(rows)}')
    return rows


def regenerate(sequence, kind, seed, cfg):
    rng = np.random.default_rng(seed)
    return dinucleotide_shuffle(sequence, rng) if kind == 'dinucleotide' else shuffled(sequence, rng, cfg['shuffle_block_nt'])


def verify_batch(task):
    sequence, kind, cfg, draws = task
    return [(seed, stored, refinement_scan(regenerate(sequence, kind, seed, cfg), cfg, False)['score']) for seed, stored in draws]


def pilercr_statistics(sequence):
    """(best-array score, best copies, crashed); a crash scores +inf, an exceedance of anything."""
    try:
        arrays = run_pilercr(sequence, PILERCR_TUNED)
    except subprocess.CalledProcessError:
        return math.inf, math.inf, True
    return best_array_score(arrays), max((a['copies'] for a in arrays), default=0), False


def pilercr_batch(task):
    sequence, kind, cfg, seeds = task
    return [(seed, *pilercr_statistics(regenerate(sequence, kind, seed, cfg))) for seed in seeds]


def run(args):
    cfg, items = search_config(), inputs()
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        # 1. The regenerated shuffles are the stored ones: D rescoring must match exactly.
        checks = [(sequence, kind, cfg, [tuple(draws[kind][i]) for i in VERIFY_DRAWS if i < len(draws[kind])])
                  for _, _, sequence, draws, _ in items for kind in KINDS]
        verified = [row for batch in pool.map(verify_batch, checks) for row in batch]
        mismatches = [row for row in verified if row[1] != row[2]]
        if mismatches:
            raise ValueError(f'{len(mismatches)} regenerated shuffles do not reproduce stored D scores: {mismatches[:3]}')
        observed_D = list(pool.map(refinement_scan, [x[2] for x in items], [cfg]*len(items)))
        for (input_id, _, _, _, stored), scan in zip(items, observed_D):
            if scan['score'] != stored:
                raise ValueError('Observed D score changed: '+input_id)
        print(f'verified {len(verified)} regenerated shuffles and {len(items)} observed D scores; '
              f'{time.perf_counter()-started:.0f} s', flush=True)
        # 2. Tuned PILER-CR on every observed window and every stored shuffle.
        observed = list(pool.map(pilercr_statistics, [x[2] for x in items]))
        labels, tasks = [], []
        for input_id, _, sequence, draws, _ in items:
            for kind in KINDS:
                seeds = [int(s) for s, _ in draws[kind]]
                for i in range(0, len(seeds), 20):
                    labels.append((input_id, kind))
                    tasks.append((sequence, kind, cfg, seeds[i:i+20]))
        nulls = {}
        for (input_id, kind), batch in zip(labels, pool.map(pilercr_batch, tasks)):
            nulls.setdefault(input_id, {}).setdefault(kind, []).extend(batch)
    print(f'PILER-CR on {sum(len(v) for d in nulls.values() for v in d.values())} shuffles; '
          f'{time.perf_counter()-started:.0f} s', flush=True)
    rows, null_rows = [], []
    for (input_id, group, sequence, draws, d_score), (t_obs, c_obs, observed_crash) in zip(items, observed):
        row = dict(input_id=input_id, group=group, window_nt=len(sequence),
                   draws_per_background=len(draws['local_blocks']), D_observed_score=d_score)
        for kind in KINDS:
            row[f'D_{kind}_p'] = empirical(d_score, [v for _, v in draws[kind]])['empirical_p']
        row['D_conservative_p'] = max(row[f'D_{k}_p'] for k in KINDS)
        if observed_crash:
            raise ValueError('PILER-CR crashed on an observed window: '+input_id)
        row.update(pilercr_array_called=t_obs > 0, pilercr_best_copies=c_obs, pilercr_observed_score=t_obs)
        for kind in KINDS:
            values = nulls[input_id][kind]
            if [s for s, _, _, _ in values] != [int(s) for s, _ in draws[kind]]:
                raise ValueError('PILER-CR null seeds are out of order: '+input_id)
            row[f'pilercr_{kind}_p'] = empirical(t_obs, [t for _, t, _, _ in values])['empirical_p']
            row[f'pilercr_{kind}_p_crash_as_no_array'] = empirical(t_obs, [0. if x else t for _, t, _, x in values])['empirical_p']
            row[f'pilercr_{kind}_crashes'] = sum(x for _, _, _, x in values)
            row[f'pilercr_{kind}_null_max'] = max((t for _, t, _, x in values if not x), default=0.)
            row[f'pilercr_{kind}_shuffles_with_array'] = sum(t > 0 and not x for _, t, _, x in values)
            row[f'pilercr_copies_{kind}_p'] = empirical(c_obs, [c for _, _, c, _ in values])['empirical_p']
            null_rows.extend(dict(input_id=input_id, background=kind, draw=i, seed=s, pilercr_score=t, pilercr_best_copies=c,
                                  pilercr_crashed=x) for i, (s, t, c, x) in enumerate(values))
        row['pilercr_conservative_p'] = max(row[f'pilercr_{k}_p'] for k in KINDS)
        row['pilercr_conservative_p_crash_as_no_array'] = max(row[f'pilercr_{k}_p_crash_as_no_array'] for k in KINDS)
        row['pilercr_copies_conservative_p'] = max(row[f'pilercr_copies_{k}_p'] for k in KINDS)
        rows.append(row)
    write_tsv(OUT/'comparison.tsv', rows, list(rows[0]))
    with gzip.open(OUT/'pilercr_null_scores.tsv.gz', 'wt') as handle:
        handle.write('\t'.join(null_rows[0])+'\n')
        handle.writelines('\t'.join(str(v) for v in r.values())+'\n' for r in null_rows)
    write_tsv(OUT/'shuffle_regeneration_checks.tsv',
              [dict(seed=s, stored_D_score=a, regenerated_D_score=b) for s, a, b in verified], ['seed', 'stored_D_score', 'regenerated_D_score'])
    m = len(rows)
    summary = dict(inputs=m, draws_per_background=rows[0]['draws_per_background'], pilercr_options=list(PILERCR_TUNED),
        statistic='PILER-CR best-array score (copies-1) x consensus repeat nt x mean copy identity; 0 if no array',
        secondary_statistic='PILER-CR best-array copy count',
        regenerated_shuffles_verified_against_stored_D_scores=len(verified),
        pilercr_crashed_shuffles=sum(r['pilercr_crashed'] for r in null_rows), elapsed_seconds=round(time.perf_counter()-started, 1))
    for group, members in (('designed_planted', [r for r in rows if r['group'] != 'published_reference']),
                           ('published_reference', [r for r in rows if r['group'] == 'published_reference']),
                           ('all', rows)):
        summary[group] = {name: dict(uncorrected_le_0_05=sum(r[name] <= .05 for r in members),
                                     bonferroni_31_le_0_05=sum(min(1, r[name]*m) <= .05 for r in members),
                                     inputs=len(members))
                          for name in ('D_conservative_p', 'pilercr_conservative_p', 'pilercr_conservative_p_crash_as_no_array',
                                       'pilercr_copies_conservative_p')}
        summary[group]['pilercr_observed_array_called'] = sum(r['pilercr_array_called'] for r in members)
    write_json(OUT/'summary.json', summary)
    print(json.dumps(summary, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=8)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
