"""Independently audit exported prediction provenance and make a portable report.

Run in Slurm after reporting. The compact archive contains tables, figures,
inputs, and provenance. All structure models remain in runs/, linked by hash.
"""
from collections import Counter
import csv
import io
import json
import os
from pathlib import Path
import shutil
import zipfile

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

from art_binding_api_credentials import credential
from art_bindings_common import ROOT,OUT,utc,sha,write,read_table,payload_hash


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    protocol=json.loads((OUT/'protocol.json').read_text())
    for name,digest in protocol['frozen_files'].items():assert sha(OUT/name)==digest,name
    for name,digest in protocol['source_sha256'].items():assert sha(ROOT/name)==digest,name
    samples={r['sample_id']:r for r in read_table(OUT/'inputs/art_samples.csv')}
    cases=read_table(OUT/'experiments.csv');index=json.loads((OUT/'request_index.json').read_text())
    assert len(cases)==protocol['cases']==1221
    assert len(index)==protocol['unique_requests']==1221
    assert dict(Counter(c['stage'] for c in cases))==protocol['case_counts_by_stage']
    requests={r['request_hash']:r for r in index};assert len(requests)==len(index)
    seen_cases=set();payloads={};fixed_inputs={};checked_input_chains=0
    for c in cases:
        assert c['case_id'] not in seen_cases;seen_cases.add(c['case_id'])
        record=requests[c['request_hash']]
        assert c['case_id'] in record['case_ids'] and c['stage'] in record['stages']
        payload=json.loads((OUT/c['request_path']).read_text())
        assert payload_hash(payload)==c['request_hash']
        payloads[c['request_hash']]=payload
        assert payload['diffusion_samples']==3 and payload['sampling_steps']==200 and payload['recycling_steps']==3
        assert payload['without_potentials'] is True and payload['concatenate_msas'] is False
        assert payload['constraints']==[] and payload['ligands']==[]
        metadata=json.loads(c['chain_metadata_json']);assert len(metadata)==len(payload['polymers'])
        for meta,polymer in zip(metadata,payload['polymers']):
            source=samples[meta['sample_id']];role=meta['role']
            original=source[role+'_protein_sequence'];start,end=meta['start0'],meta['end0']
            assert 0<=start<end<=len(original)
            if role=='RT':
                b=int(source['NTD_length_aa']);assert meta['source_NTD_end0']==b
                expected={'full':(0,len(original)),'NTD_baseline':(0,b),'NTD_minus20':(0,b-20),
                    'NTD_plus20':(0,b+20),'core_only':(b,len(original))}[c['representation']]
                assert (start,end)==expected
                assert meta['sample_id']==c['RT_sample_id']
            else:
                assert (start,end)==(0,len(original))
                assert meta['sample_id']==c['partner_sample_id']
            assert polymer['id']==meta['chain_id'] and polymer['sequence']==original[start:end]
            assert polymer['molecule_type']=='protein'
            msa=polymer['msa']['controlled_input']
            if c['condition']=='S_query_only':
                assert msa['a3m']['alignment']=='>query\n'+polymer['sequence']+'\n'
            else:
                alignment=msa['csv']['alignment'];reader=csv.DictReader(io.StringIO(alignment))
                assert reader.fieldnames==['key','sequence'];rows=list(reader)
                assert rows[0]['sequence']==polymer['sequence'] and not any(r['key'] for r in rows)
                assert all(len(''.join(x for x in r['sequence'] if not x.islower()))==len(polymer['sequence']) for r in rows)
            # A chain's input must not change with its pairing or source alias.
            fixed_key=(source['partner_type'],role,polymer['sequence'],c['condition'])
            assert fixed_inputs.setdefault(fixed_key,msa)==msa
            checked_input_chains+=1
        if len(metadata)==2:
            assert (c['native_pair']=='True')==(c['RT_sample_id']==c['partner_sample_id'])
        else:assert c['native_pair']=='False'
        assert c['confirmatory']=='False'
    assert {x for r in index for x in r['case_ids']}==seen_cases
    # Validate the immutable reporting snapshot while the inference launcher
    # may continue completing additional requests in the background.
    scores=read_table(OUT/'results/model_scores.csv');medians=read_table(OUT/'results/pair_medians.csv')
    exported_hashes={r['request_hash'] for r in medians}
    states=Counter();models=0;chain_checks=0;independent_sequence_checks=0;complete_hashes=set()
    for digest,record in requests.items():
        folder=OUT/'runs'/digest;state_path=folder/'state.json'
        state=json.loads(state_path.read_text()) if state_path.exists() else {'status':'not_started'}
        states[state['status']]+=1
        if digest not in exported_hashes:continue
        assert state['status']=='complete'
        normalized=json.loads((folder/'normalized.json').read_text());payload=payloads[digest]
        assert normalized['verified'] and len(normalized['models'])==3
        assert normalized['verified_sequences']==3*len(payload['polymers'])
        for model in normalized['models']:
            assert sha(folder/model['structure_file'])==model['structure_sha256']
            assert (folder/model['coordinates_file']).is_file()
            if len(payload['polymers'])==2:assert (folder/model['contacts_file']).is_file()
            models+=1
        # Replay the sequence check using a different parser and the mmCIF
        # label identifiers, independently of normalize_response's PDB objects.
        raw=MMCIF2Dict(str(folder/normalized['models'][0]['structure_file']))
        residues={}
        for chain,residue,resname,atom in zip(raw['_atom_site.label_asym_id'],raw['_atom_site.label_seq_id'],
                raw['_atom_site.label_comp_id'],raw['_atom_site.label_atom_id']):
            if atom=='CA':residues.setdefault(chain,{})[int(residue)]=resname
        assert set(residues)=={p['id'] for p in payload['polymers']}
        for p in payload['polymers']:
            sequence=''.join(seq1(name) for i,name in sorted(residues[p['id']].items()))
            assert sequence==p['sequence'];independent_sequence_checks+=1
        chain_checks+=normalized['verified_sequences'];complete_hashes.add(digest)
    assert len(scores)==models and len(medians)==len(complete_hashes)
    for row in scores:
        assert row['request_hash'] in complete_hashes and sha(OUT/row['structure_file'])==row['structure_sha256']
        assert row['confirmatory']=='False'
    for row in read_table(OUT/'results/matrix_contrasts.csv'):
        if row['complete_matrix']!='True':assert row['descriptive_diagonal_minus_offdiagonal']==''
        assert row['confirmatory']=='False' and row['pvalue']==''
    complete=len(complete_hashes)==len(index)
    result=dict(status='passed' if complete else 'passed_partial',updated_utc=utc(),complete_prediction_batch=complete,
        cases_checked=len(cases),unique_requests_checked=len(index),input_chain_source_checks=checked_input_chains,
        completed_requests=len(complete_hashes),completed_models=models,normalized_model_chain_checks=chain_checks,
        independent_mmCIF_chain_checks=independent_sequence_checks,request_statuses=dict(states),
        biological_binding_validated=False,confirmatory=False,slurm_job_id=os.environ['SLURM_JOB_ID'])
    shutil.copy2(Path(__file__),OUT/'provenance/verification_code.py')
    write(OUT/'results/verification.json',result)
    write(OUT/'checkpoints/05_results_verified.json',result)
    files=[OUT/name for name in ['README.md','protocol.json','experiments.csv','request_index.json']]
    files.extend(p for directory in ['inputs','panels','results','figures','checkpoints','provenance']
        for p in (OUT/directory).rglob('*') if p.is_file())
    key,provider,mode=credential(ROOT)
    for path in files:
        assert key.encode() not in path.read_bytes(),'Credential found in an export; archive not written'
    result['export_credential_scan_passed']=True
    write(OUT/'results/verification.json',result);write(OUT/'checkpoints/05_results_verified.json',result)
    hashes={str(p.relative_to(OUT)):sha(p) for p in files}
    write(OUT/'paper_results_manifest.json',hashes)
    archive=OUT/'paper_results.zip';temporary=archive.with_suffix('.zip.tmp')
    with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path in files+[OUT/'paper_results_manifest.json']:
            z.write(path,'results_bindings/'+str(path.relative_to(OUT)))
        z.writestr('results_bindings/ARCHIVE_SCOPE.txt',
            'This compact archive contains inputs, tables, figures and provenance.\n'
            'Full mmCIF models, per-residue contacts and raw API responses remain under runs/ in the source project.\n'
            'Model paths and immutable SHA256 hashes are recorded in results/model_scores.csv.\n'
            'Check results/verification.json for partial versus complete coverage. These are computational predictions.\n')
    temporary.replace(archive)
    write(OUT/'paper_results_archive.json',dict(file=archive.name,sha256=sha(archive),bytes=archive.stat().st_size,
        complete_prediction_batch=complete,created_utc=utc()))
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
