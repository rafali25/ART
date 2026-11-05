"""Close the bounded five-step run without changing its scientific checkpoints."""
from collections import defaultdict
from importlib.metadata import version
import hashlib
import json
import os
from pathlib import Path
import shutil
import socket
import sys

from Bio import SeqIO
from Bio.Seq import Seq

from art_expansion_workflow import ROOT, OUT, PROFILES, DATASETS, sha, file_md5, read_tsv, table, write_json, utc


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Use a Slurm CPU batch job; login-node analysis is prohibited.')
    audit = OUT/'audit'
    audit.mkdir(exist_ok=True)
    checkpoints = [json.loads((OUT/f'checkpoints/{n:02d}.json').read_text()) for n in range(1, 6)]
    checks = []
    for n, checkpoint in enumerate(checkpoints, 1):
        assert checkpoint['step'] == n and checkpoint['status'] == 'complete'
        for name, expected in checkpoint['artifacts'].items():
            assert sha(ROOT/name) == expected, name
            checks.append(dict(scope=f'checkpoint_{n}', path=name, sha256=expected, passed=True))
    baseline = ROOT/'results/art_type1_catalog'
    baseline_summary = json.loads((baseline/'summary.json').read_text())
    for base, manifest in [(baseline, json.loads((baseline/'output_sha256.json').read_text())),
                           (ROOT, baseline_summary['input_manifest_sha256'])]:
        for name, expected in manifest.items():
            path = base/name.replace('\\', '/')
            assert sha(path) == expected, str(path)
            checks.append(dict(scope='frozen_baseline', path=str(path.relative_to(ROOT)), sha256=expected, passed=True))
    policy = json.loads((OUT/'cluster_execution.json').read_text())
    assert sha(ROOT/policy['policy_file']) == policy['policy_sha256']
    downloads = []
    for name, metadata in DATASETS.items():
        path = OUT/'downloads'/metadata['file']
        assert path.stat().st_size == metadata['bytes']
        assert file_md5(path) == metadata['md5']
        downloads.append(dict(dataset=name, bytes=path.stat().st_size, md5=metadata['md5'], sha256=sha(path)))
    table(audit/'verified_artifacts.tsv', checks)
    table(audit/'download_integrity.tsv', downloads)

    # The initial deduplication registry omitted the historical type-II source
    # manifest. Audit that omission explicitly; final locus dedup already removes
    # gene-model aliases, so these screen counters must not be called novelty.
    old_sources = defaultdict(list)
    for row in read_tsv(ROOT/'results/art_expansion/stage1/source_records.tsv'):
        path = ROOT/row['path'].replace('\\', '/')
        assert sha(path) == row['file_sha256']
        record = SeqIO.read(path, row['format'])
        dna = str(record.seq).upper()
        canonical = min(hashlib.sha256(dna.encode()).hexdigest(),
                        hashlib.sha256(str(Seq(dna).reverse_complement()).encode()).hexdigest())
        old_sources[canonical].append(row['accession'])
    overlaps = []
    for dataset in ['HPGC', 'JPGD', 'UHGV_selected']:
        for row in read_tsv(OUT/dataset/'genome_manifest.tsv'):
            for accession in old_sources.get(row.get('canonical_sequence_sha256'), []):
                overlaps.append(dict(dataset=dataset, genome_key=row['genome_key'],
                    original_id=row['original_id'], screen_status=row['status'],
                    prior_source_accession=accession,
                    canonical_sequence_sha256=row['canonical_sequence_sha256']))
    table(audit/'prior_source_overlap.tsv', overlaps)

    rows = read_tsv(OUT/'step5/pair_catalog.tsv')
    aliases = read_tsv(OUT/'step5/exact_pair_aliases.tsv')
    verification = json.loads((OUT/'step5/verification.json').read_text())
    summary = json.loads((OUT/'step5/summary.json').read_text())
    assert verification['passed'] and verification['pair_source_maps_replayed'] == 49
    assert len(rows) == summary['analysis_locus_representatives'] == 46
    assert len(aliases) == summary['pair_or_gene_model_aliases'] == 3
    assert len({r['pair_sequence_sha256'] for r in rows}) == len(rows)
    cohort = []
    for dataset in sorted({r['dataset'] for r in rows}):
        selected = [r for r in rows if r['dataset'] == dataset]
        type_i = [r for r in selected if r['family'] == 'Type_I']
        complete = [r for r in type_i if r['pair_complete'] == 'True']
        cohort.append(dict(dataset=dataset, typeI_pairs=len(type_i), complete_typeI_pairs=len(complete),
            complete_strict_typeI_pairs=sum(r['ART_profile_passes'] == 'True' for r in complete),
            typeII_pairs=sum(r['family'] == 'Type_II' for r in selected)))
    table(audit/'cohort_counts.tsv', cohort)
    strict = read_tsv(OUT/'step5/complete_typeI_strict_screen_pairs.tsv')
    complete = read_tsv(OUT/'step5/complete_typeI_pairs.tsv')
    assert len(strict) == summary['complete_strict_screen_typeI_pairs'] == 18
    assert len(complete) == summary['complete_typeI_pairs'] == 31
    edge_rows = [r for r in read_tsv(OUT/'step5/matching_pair_comparisons.tsv')
                 if r['passes_sequence_matching'] == 'True']
    assert len(edge_rows) == summary['block_feasibility']['matching_pair_edges'] == 1
    assert len(read_tsv(OUT/'step5/all_matching_triples.tsv')) == summary['block_feasibility']['matching_triples'] == 0
    assert summary['block_feasibility']['confirmatory_structural_design_ready'] is False

    # Preserve the exact closing code and source helpers beside the evidence.
    code = OUT/'code'
    code.mkdir(exist_ok=True)
    paths = sorted((ROOT/'scripts').glob('art_expansion_*'))
    paths += [ROOT/'scripts/test_art_expansion_workflow.py', ROOT/'scripts/test_art_type1_public.py',
              ROOT/'scripts/art_type1_pipeline.py', ROOT/'scripts/art_reconstruction_common.py',
              ROOT/'scripts/extract_art_loci.py', ROOT/'scripts/art_type1_public.py']
    for path in paths:
        if path.is_file():
            shutil.copy2(path, code/path.name)
    native = ROOT/'results/art_sample_expansion_20261104/deferred_acquisition/screen_deferred.py'
    shutil.copy2(native, code/'prior_screen_deferred.py')
    write_json(audit/'software.json', dict(python=sys.version, packages={
        p: version(p) for p in ['biopython', 'numpy', 'pyrodigal-gv', 'zstandard']},
        profiles={name: sha(path) for name, path in PROFILES.items()},
        source_scripts={str(path.relative_to(ROOT)): sha(path) for path in paths if path.is_file()},
        native_wrapper_sha256=sha(native)))
    logs = OUT/'logs'
    logs.mkdir(exist_ok=True)
    log_names = ['art-prepare-371266.log', 'art-download-371284.log', 'art-leads-371364.log',
                 'art-catalogs-371321.log', 'art-members-371425.log',
                 'art-feasibility-371381.log', 'art-feasibility-371557.log']
    copied_logs = []
    for name in log_names:
        path = ROOT/'tmp'/name
        if path.exists():
            shutil.copy2(path, logs/name)
            copied_logs.append(name)
    test_log = (ROOT/'tmp/art-feasibility-371557.log').read_text()
    assert 'Ran 4 tests' in test_log and 'Ran 8 tests' in test_log and test_log.count('\nOK\n') == 2

    (OUT/'README.md').write_text('''# Completed five-step ART acquisition run

This bounded run adds **21 distinct computational type I candidate pairs**, of
which **19 have complete CDSs**, to the pre-run collection. The merged collection
has **36 type I pairs, 31 complete pairs, and 18 complete pairs passing the frozen
strict ART screen**. It also retains 10 type II pairs as a separate arm.
These counts do not establish biological function or independent samples.

| Step | Executed work | Checkpoint |
|---|---|---|
| 1 | Reviewed four prior additions; retained strict and review tiers; verified 24 frozen baseline hashes | [01](checkpoints/01.md) |
| 2 | Attempted all 45 leads from 44 runs; exactly mapped 24 source contigs; screened them and queried MGnify occurrences/members | [02](checkpoints/02.md) |
| 3 | Downloaded and checksum-verified all 7,230 HPGC records; called/searched 7,151 records after initial DNA deduplication | [03](checkpoints/03.md) |
| 4 | Downloaded and checksum-verified all 10,754 JPGD records; called/searched 4,853 after initial deduplication; retrieved 107 selected UHGV genomes and screened 66 additional DNA sequences | [04](checkpoints/04.md) |
| 5 | Replayed source maps and CDS translations; built separate RT-core/partner trees; deduplicated pairs; evaluated matched blocks | [05](checkpoints/05.md) |

HPGC contributed one distinct complete type I pair. Its type II hit is the old
Machias locus with a different partner gene boundary. JPGD contributed four
distinct complete type I pairs, including one strict ART-screen pair; two other
JPGD pair records were exact duplicates. The archive screen contributed 16
distinct type I pairs. Selected UHGV screening found no ART/partner pair.
See [cohort counts](audit/cohort_counts.tsv) and [aliases](step5/exact_pair_aliases.tsv).

The 36 type I pairs span 22 partner clusters at 70% identity/80% mutual coverage;
the complete-pair subset has 19 such representatives. Cluster counts are distinct
from complete strict-screen samples and from phylogenetically independent loci.

The strict-screen matching graph has **one compatible pair and zero matching
triples**, hence **zero disjoint three-locus blocks** under the proposed criteria
(50–90% identity in both RT core and partner; >=80% mutual coverage; compatible
lengths; distinct source groups). Allowing tree-affinity review candidates gives
six compatible pairs but still zero triples. The proposed eight-block minimum
is unmet. Architecture, phylogenetic independence, calibration, power, and RNA
eligibility remain additional gates. See [matching comparisons](step5/matching_pair_comparisons.tsv)
and [summary](step5/summary.json).

The RT-core tree contains 137 sequences including references; the partner tree
contains 36. FastTree internal values are local SH-like support, not bootstraps.
Tree affinity does not validate low-coverage ART hits or replace missing NART
controls. No new array-significance tests or accepted arrays were produced.

All computation, downloads, and tests ran in Slurm CPU batch jobs under
`cscc-cpu-p`/`cscc-cpu-qos`, following the supplied cluster document. Twelve
mapping tests passed; all 49 raw pair source maps and 98 CDS translations were
replayed successfully. [Final verification](verification.json),
[artifact checks](audit/verified_artifacts.tsv), and [code snapshot](code/) retain
the closing evidence. Original catalogs, profiles, nulls, and significance
results were preserved.

Retrieval limitations remain explicit: 21 archive leads were unresolved within
the 200-MB successful-transfer cap; the MGnify remote Parquet endpoint refused
connections (zero data bytes), so member expansion remains unresolved; 13 of 120
selected UHGV genomes could not be retrieved. These are not absence calls.
A separate IMG/VR download was not performed. The initial DNA registry omitted
the historical type-II source manifest; [source-overlap audit](audit/prior_source_overlap.tsv)
records that overlap, and final locus deduplication removes the Machias alias.
The 7,151/4,853 counts describe records called in the initial registry, not novel
genomes. Slurm accounting (`sacct`) was unavailable because its database
connection was refused; scheduler output and compute logs retain run evidence.

The final SHA-256 manifest covers every durable file in this run except itself,
`verification.json`, and temporary files. Checkpoint revisions are retained in
`checkpoints/history/`; each step's recorded artifact hashes were rechecked.
''')
    with (OUT/'RUN_NOTES.md').open('a') as handle:
        handle.write('''

Final closure:

- 371381 passed all 12 tests and source-map replay but stopped because MMseqs
  rejected the spelling `.7`/`.8` for numeric thresholds. Arguments were changed
  to `0.7`/`0.8`, with the same biological thresholds.
- 371557 completed the rerun and saved checkpoint 5. No source acquisition or
  bulk gene searches were repeated for this correction.
- The bounded Parquet member attempt in 371425 failed on an official endpoint
  connection refusal, transferring zero Parquet bytes. No member absence is
  inferred. Source-run ENA metadata retrieval and archive screening completed.
- Final audit runs as a separate CPU Slurm job, preserves a source snapshot,
  rechecks all checkpoint artifacts and original baseline hashes, and writes
  `verification.json`. No active analysis jobs remain after closure.
- `sacct` could not reach its accounting database. This is recorded rather than
  substituting invented accounting statistics.
''')
    result = dict(passed=True, audited_utc=utc(), slurm_job_id=os.environ['SLURM_JOB_ID'],
        compute_host=socket.gethostname(), checkpoints_complete=5,
        checkpoint_artifact_hashes_verified=sum(r['scope'].startswith('checkpoint') for r in checks),
        original_catalog_hashes_unchanged=sum(r['scope'] == 'frozen_baseline' for r in checks),
        published_catalog_MD5s_verified=2, policy_sha256_verified=True, tests_passed=12,
        pair_source_maps_replayed=49, CDS_translations_replayed=98,
        final_pair_records=46, pair_or_gene_model_aliases=3, typeI_pairs=36,
        complete_typeI_pairs=31, complete_strict_typeI_pairs=18,
        current_run_distinct_typeI_pairs=21, current_run_complete_typeI_pairs=19,
        disjoint_sequence_matched_blocks=0, confirmatory_structural_design_ready=False,
        historical_typeII_source_overlap_records=len(overlaps),
        logs_preserved=copied_logs, unresolved_archive_leads=21,
        parquet_member_retrieval='connection_refused_zero_bytes_not_absence',
        selected_UHGV_retrieval_failures=13, slurm_accounting='database_connection_refused',
        artifact_manifest='output_sha256.json',
        manifest_excludes=['output_sha256.json', 'verification.json', '*.tmp'])
    manifest = {str(path.relative_to(OUT)): sha(path) for path in sorted(OUT.rglob('*'))
                if path.is_file() and path not in [OUT/'output_sha256.json', OUT/'verification.json']
                and not path.name.endswith('.tmp')}
    write_json(OUT/'output_sha256.json', manifest)
    for name, expected in manifest.items():
        assert sha(OUT/name) == expected, name
    result['durable_artifact_hashes_verified'] = len(manifest)
    write_json(OUT/'verification.json', result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
