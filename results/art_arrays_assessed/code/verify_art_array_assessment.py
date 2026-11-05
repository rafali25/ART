"""Check cohort correction, seeded null replay and original-source DNA mappings."""
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

from art_array_assessment import ROOT,dinucleotide_shuffle,sensitivity_batch
from art_array_detector import shuffled
from extract_art_loci import read_tsv,write_json


def verify():
    output = ROOT/'results/art_arrays_assessed'
    inputs = ROOT/'results/oriented_art_loci'
    config = json.loads((output/'run_config.json').read_text())
    cfg,assessment = config['search'],config['assessment']
    manifests = {r['locus_id']:r for r in read_tsv(inputs/'locus_manifest.tsv')}
    sources = {r['accession']:r for r in read_tsv(ROOT/'results/art_discovery/source_records.tsv')}
    originals = {}
    for accession in {r['accession'] for r in manifests.values()}:
        source = sources[accession]
        path = ROOT/source['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == source['file_sha256']
        with path.open(encoding='utf-8') as handle:
            originals[accession] = next(r for r in SeqIO.parse(handle,source['format']) if r.id == accession)
    loci = {l:SeqIO.read(inputs/r['genbank_path'],'genbank') for l,r in manifests.items()}
    groups = defaultdict(list)
    for row in read_tsv(output/'null_best_scores.tsv'):
        groups[row['locus_id'],row['background'],row['window_nt']].append(row)
    replayed = 0
    for (locus,kind,window),rows in groups.items():
        expected_n = assessment['dinucleotide_shuffles_per_window' if kind == 'dinucleotide' else 'local_inference_shuffles_per_window']
        assert len(rows) == expected_n
        assert sorted(int(r['shuffle_index']) for r in rows) == list(range(expected_n))
        rt = int(manifests[locus]['RT_locus_start0'])
        sequence = str(loci[locus].seq[max(0,rt-int(window)):rt]).upper()
        sample = rows[len(rows)//2]
        seed = int(sample['rng_seed'])
        rng = np.random.default_rng(seed)
        permuted = dinucleotide_shuffle(sequence,rng) if kind == 'dinucleotide' else shuffled(sequence,rng,cfg['shuffle_block_nt'])
        if kind == 'dinucleotide':
            assert Counter(zip(sequence,sequence[1:])) == Counter(zip(permuted,permuted[1:]))
            assert Counter(sequence) == Counter(permuted)
        else:
            for start in range(0,len(sequence),cfg['shuffle_block_nt']):
                assert Counter(sequence[start:start+50]) == Counter(permuted[start:start+50])
        actual = sensitivity_batch((sequence,kind,cfg,[seed]))[0][1]
        assert abs(actual-float(sample['best_chain_score'])) < 1e-8
        replayed += 1
    background = read_tsv(output/'background_tests.tsv')
    assert len(background) == len(manifests)
    for row in background:
        locus = row['locus_id']
        for kind in ('local_blocks','dinucleotide'):
            for window in ('3000','6000'):
                observed = float(row[kind+'_'+window+'_best_score'])
                values = [float(r['best_chain_score']) for r in groups[locus,kind,window]]
                p = (1+sum(s >= observed for s in values))/(1+len(values))
                assert abs(p-float(row[kind+'_'+window+'_p'])) < 1e-12
        primary = max(float(row['local_blocks_6000_p']),float(row['dinucleotide_6000_p']))
        assert float(row['bonferroni_p']) == min(1.,primary*len(manifests))
    array_rows = read_tsv(output/'array_assessments.tsv')
    for row in array_rows:
        locus = row['locus_id']
        observed = float(row['seed_chain_score'])
        primary_p = []
        for kind in ('local_blocks','dinucleotide'):
            values = [float(r['best_chain_score']) for r in groups[locus,kind,'6000']]
            p = (1+sum(s >= observed for s in values))/(1+len(values))
            assert abs(p-float(row[kind+'_6000_p'])) < 1e-12
            primary_p.append(p)
        assert float(row['bonferroni_p']) == min(1.,max(primary_p)*len(manifests))
        if row['assessment_status'] == 'accepted_array_candidate':
            assert int(row['copy_count']) >= 3
            assert float(row['bonferroni_p']) <= assessment['familywise_alpha']
            assert not row['coding_overlap_ids'] and not row['intervening_gene_ids']
            assert int(row['ambiguous_array_nt']) == 0
        details = json.loads((output/'details'/f"{row['array_id']}.json").read_text())
        training = output/'alignments'/f"{row['array_id']}.training"
        with Path(str(training)+'.flanks.fna').open() as handle:
            raw = {r.id:str(r.seq).upper() for r in SeqIO.parse(handle,'fasta')}
        with Path(str(training)+'.mafft.fna').open() as handle:
            aligned = {r.id:str(r.seq).upper() for r in SeqIO.parse(handle,'fasta')}
        assert raw == {key:text.replace('-','') for key,text in aligned.items()}
        rt = int(manifests[locus]['RT_locus_start0'])
        seq = str(loci[locus].seq[max(0,rt-cfg['extended_window_nt']):rt]).upper()
        for i,p in enumerate(details['chain']['positions'],1):
            assert raw[f'copy{i}'] == seq[max(0,p-cfg['boundary_flank_nt']):min(len(seq),p+cfg['refinement_seed_nt']+cfg['boundary_flank_nt'])]
    def original_text(locus,pieces):
        source = originals[manifests[locus]['accession']]
        parts = []
        for p in pieces:
            text = source.seq[p['source_start0']:p['source_end0']]
            parts.append(text if p['orientation'] == 1 else text.reverse_complement())
        return str(Seq('').join(parts))
    sequence_count = 0
    all_fasta = {}
    for table,fasta,key in (('repeats.tsv','repeats.fna','repeat_id'),('spacers.tsv','spacers.fna','spacer_id')):
        with (output/fasta).open() as handle:
            records = {r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
        all_fasta.update(records)
        for row in read_tsv(output/table):
            locus = row['locus_id']
            a,b = int(row['locus_start0']),int(row['locus_end0'])
            assert records[row[key]] == row['sequence'] == str(loci[locus].seq[a:b])
            assert original_text(locus,json.loads(row['source_segments_json'])) == row['sequence']
            if table == 'spacers.tsv':
                assert int(row['length_nt']) == int(row['repeat_start_spacing_nt'])-int(row['preceding_repeat_length_nt'])
            sequence_count += 1
    with (output/'arrays.fna').open() as handle:
        array_sequences = {r.id:str(r.seq) for r in SeqIO.parse(handle,'fasta')}
    all_fasta.update(array_sequences)
    for row in array_rows:
        if row['array_id'] not in array_sequences: continue
        expected = original_text(row['locus_id'],json.loads(row['source_segments_json']))
        assert array_sequences[row['array_id']] == expected
        sequence_count += 1
    accepted = {r['array_id'] for r in array_rows if r['assessment_status'] == 'accepted_array_candidate'}
    for name in ('accepted_arrays.fna','accepted_repeats.fna','accepted_spacers.fna'):
        with (output/name).open() as handle:
            for record in SeqIO.parse(handle,'fasta'):
                assert str(record.seq) == all_fasta[record.id]
                array_id = record.id.rsplit('.repeat',1)[0].rsplit('.spacer',1)[0]
                assert array_id in accepted
    exploratory = read_tsv(output/'exploratory_two_copy.tsv')
    for row in exploratory:
        assert row['assessment_status'] == 'exploratory_two_copy_not_accepted'
        assert row['copy_count'] == '2' and not row['calibrated_p'] and not row['bonferroni_p']
        assert row['array_id'] not in accepted
        for text,pieces in zip(json.loads(row['repeat_sequences_json']),json.loads(row['source_segments_json'])):
            assert original_text(row['locus_id'],pieces) == text
            sequence_count += 1
    result = dict(passed=True,loci_in_testing_family=len(manifests),null_draws_checked=sum(map(len,groups.values())),
        independent_null_groups_replayed=replayed,sequence_mappings_verified=sequence_count,
        repeat_training_alignments_verified=len(array_rows),two_copy_labels_verified=len(exploratory),
        accepted_arrays_verified=len(accepted),all_original_sources_unchanged=True,
        all_empirical_p_values_and_cohort_corrections_verified=True)
    write_json(output/'verification.json',result)
    shutil.copyfile(__file__,output/'code'/Path(__file__).name)
    write_json(output/'output_sha256.json',{p.relative_to(output).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'output_sha256.json'})
    for relative,digest in json.loads((output/'output_sha256.json').read_text()).items():
        assert hashlib.sha256((output/relative).read_bytes()).hexdigest() == digest
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    verify()
