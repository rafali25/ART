"""Independent output/source checks and sampled replay of seeded null scans."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

from art_array_detector import ROOT, null_batch, pwm_scores, shuffled
from extract_art_loci import read_tsv, write_json


def verify(output=ROOT/'results/art_arrays',inputs=ROOT/'results/oriented_art_loci'):
    provenance = json.loads((output/'provenance.json').read_text())
    for path,digest in provenance['input_hashes'].items():
        assert hashlib.sha256((inputs/path).read_bytes()).hexdigest() == digest,path
    manifests = {r['locus_id']:r for r in read_tsv(inputs/'locus_manifest.tsv')}
    sources = {r['accession']:r for r in read_tsv(ROOT/'results/art_discovery/source_records.tsv')}
    original = {}
    for accession in {m['accession'] for m in manifests.values()}:
        source = sources[accession]
        path = ROOT/source['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source['file_sha256']
        with path.open(encoding='utf-8') as handle:
            original[accession] = next(r for r in SeqIO.parse(handle,source['format']) if r.id == accession)
    oriented = {l:SeqIO.read(inputs/m['genbank_path'],'genbank') for l,m in manifests.items()}
    checks = defaultdict(int)
    for table,fasta in (('repeats.tsv','repeats.fna'),('spacers.tsv','spacers.fna'),
                        ('initial_repeats.tsv','initial_repeats.fna'),('initial_spacers.tsv','initial_spacers.fna')):
        with (output/fasta).open() as handle:
            seqs = {r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
        for row in read_tsv(output/table):
            l = row['locus_id']
            a,b = int(row['locus_start0']),int(row['locus_end0'])
            identifier = row.get('repeat_id',row.get('spacer_id'))
            assert str(oriented[l].seq[a:b]) == row['sequence'] == seqs[identifier]
            source = original[manifests[l]['accession']]
            parts = []
            for piece in json.loads(row['source_segments_json']):
                part = source.seq[piece['source_start0']:piece['source_end0']]
                parts.append(part if piece['orientation'] == 1 else part.reverse_complement())
            assert str(Seq('').join(parts)) == row['sequence'],identifier
            if 'spacer_index' in row:
                assert b-a == int(row['repeat_start_spacing_nt'])-int(row['preceding_repeat_length_nt'])
            checks[table] += 1
    with (output/'arrays.fna').open() as handle:
        array_seqs = {r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    with (output/'initial_arrays.fna').open() as handle:
        array_seqs.update({r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')})
    for row in read_tsv(output/'array_sequence_map.tsv'):
        l = row['locus_id']
        a,b = int(row['locus_start0']),int(row['locus_end0'])
        expected = str(oriented[l].seq[a:b])
        assert expected == array_seqs[row['array_id']]
        source = original[manifests[l]['accession']]
        parts = []
        for p in json.loads(row['source_segments_json']):
            part = source.seq[p['source_start0']:p['source_end0']]
            parts.append(part if p['orientation'] == 1 else part.reverse_complement())
        assert str(Seq('').join(parts)) == expected
        checks['complete_array_tracts'] += 1
    grouped = defaultdict(list)
    for row in read_tsv(output/'null_scores.tsv'):
        grouped[row['locus_id'],row['stage']].append(row)
    cfg = json.loads((output/'run_config.json').read_text())['art_array_search']
    replayed = 0
    for (l,stage),rows in grouped.items():
        rt_start = int(manifests[l]['RT_locus_start0'])
        width = int(rows[0]['window_nt'])
        sequence = str(oriented[l].seq[rt_start-width:rt_start]).upper()
        # Replay one independently seeded draw from every null group.
        row = rows[len(rows)//2]
        seed,expected = int(row['rng_seed']),float(row['score'])
        if stage.startswith('PWM_'):
            array_id = stage.removeprefix('PWM_')
            details = json.loads((output/'details'/f'{array_id}.json').read_text())
            matrix = np.array(details['pwm']['pwm'])
            scores = pwm_scores(shuffled(sequence,np.random.default_rng(seed),cfg['shuffle_block_nt']),matrix)
            actual = float(scores.max())
            assert len(rows) == cfg['pwm_shuffles']
        else:
            mode = 'initial' if stage.startswith('initial_') else 'refinement'
            actual = null_batch((sequence,mode,cfg,[seed]))[0][1]
        assert abs(actual-expected) < 1e-8,(l,stage,actual,expected)
        replayed += 1
    for call in read_tsv(output/'locus_calls.tsv'):
        l = call['locus_id']
        assert call['extended_search_executed'] == 'True'
        assert len(grouped[l,'initial_mononucleotide']) == cfg['initial_shuffles']
        for label in ('3000','6000'):
            rows = grouped[l,'refinement_local_blocks_'+label]
            assert len(rows) == int(call['refinement_shuffles_per_window'])
            assert max(float(r['score']) for r in rows) == float(call['refinement_null_'+label+'_max'])
        newer = SeqIO.read(output/'loci'/f'{l}.arrays.gb','genbank')
        assert newer.seq == oriented[l].seq
        assert [(str(f.location),f.qualifiers) for f in newer.features[:len(oriented[l].features)]] == [
                (str(f.location),f.qualifiers) for f in oriented[l].features]
    for array in read_tsv(output/'arrays.tsv'):
        l = array['locus_id']
        assert int(array['copy_count']) >= cfg['minimum_copies_for_main_call']
        assert float(array['chain_score']) > max(float(array['null_3000_max']),float(array['null_6000_max']))
        lengths = json.loads(array['repeat_lengths_nt_json'])
        gaps = json.loads(array['spacer_lengths_nt_json'])
        starts = json.loads(array['repeat_start_spacings_nt_json'])
        assert gaps == [d-length for d,length in zip(starts,lengths)]
    result = dict(passed=True,sequence_rows_verified=dict(checks),
        independent_null_groups_replayed=replayed,all_original_contig_sequences_verified=True,
        all_original_annotations_retained=True,all_reported_null_maxima_verified=True,
        all_loci_extended_searched=True,spacer_formula_verified=True)
    write_json(output/'independent_verification.json',result)
    for script in (Path(__file__),Path(__file__).with_name('test_art_array_detector.py')):
        shutil.copyfile(script,output/'code'/script.name)
    write_json(output/'output_sha256.json',{p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'output_sha256.json'})
    for path,digest in json.loads((output/'output_sha256.json').read_text()).items():
        assert hashlib.sha256((output/path).read_bytes()).hexdigest() == digest
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    verify()
