"""Freeze expanded exploratory matrices, domain controls, and API payloads."""
from collections import Counter
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import random
import shutil
import subprocess

from Bio import AlignIO,SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from art_bindings_common import ROOT,OUT,utc,sha,write,table,read_table,cropped_a3m,msa_csv,polymer,payload_hash


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    source_config=ROOT/'config/art_bindings_exploratory_20261105.json'
    config=json.loads(source_config.read_text())
    if (OUT/'protocol.json').exists():
        old=json.loads((OUT/'protocol.json').read_text())
        assert old['source_config_sha256']==sha(source_config),'Frozen design changed; choose a new run version'
        print('Existing frozen design retained',flush=True);return
    for directory in ['inputs','panels','alignments','msas','requests','runs','results','figures','checkpoints','provenance/code','provenance/prior_pilots','logs']:
        (OUT/directory).mkdir(parents=True,exist_ok=True)
    source=ROOT/config['source_samples'];samples=read_table(source)
    strict=[r for r in samples if r['pair_complete']=='True' and r['category1_ART_profile_pass']=='True']
    typeI=sorted([r for r in strict if r['partner_type']=='Type_I'],key=lambda r:r['sample_id'])
    typeII=sorted([r for r in strict if r['partner_type']=='Type_II'],key=lambda r:r['sample_id'])
    review=sorted([r for r in samples if r['partner_type']=='Type_I' and r['pair_complete']=='True' and r['category1_ART_profile_pass']!='True'],key=lambda r:r['sample_id'])
    excluded=[dict(sample_id=r['sample_id'],locus_id=r['locus_id'],reason='incomplete_or_uncertain_CDS') for r in samples if r['pair_complete']!='True']
    assert (len(samples),len(typeI),len(typeII),len(review),len(excluded))==(46,18,10,13,5)
    shutil.copy2(source,OUT/'inputs/art_samples.csv');table(OUT/'panels/excluded_samples.csv',excluded)
    for name,rows in [('typeI_strict',typeI),('typeII_complete',typeII),('typeI_family_review',review)]:table(OUT/'panels'/(name+'.csv'),rows)
    msa_lookup={};canonical_msa={}
    for family,group in [('Type_I',typeI),('Type_II',typeII)]:
        for role,column in [('RT','RT_protein_sequence'),('partner','partner_protein_sequence')]:
            raw=OUT/'alignments'/(family+'_'+role+'.faa');aligned=raw.with_suffix('.aligned.faa')
            SeqIO.write([SeqRecord(Seq(r[column]),id=r['sample_id'],description='') for r in group],raw,'fasta')
            with aligned.open('w') as stdout,(OUT/'logs'/(family+'_'+role+'.mafft.log')).open('w') as stderr:
                subprocess.run([str(ROOT/'.venv/bin/mafft'),'--localpair','--maxiterate','1000','--thread','4','--inputorder',str(raw)],stdout=stdout,stderr=stderr,check=True)
            alignment={r.id:str(r.seq) for r in AlignIO.read(aligned,'fasta')}
            for query in group:
                sid=query['sample_id'];intervals={'full':(0,len(query[column]))}
                if family=='Type_I' and role=='RT':
                    b=int(query['NTD_length_aa']);n=len(query[column])
                    intervals.update(NTD_baseline=(0,b),NTD_minus20=(0,b-20),NTD_plus20=(0,b+20),core_only=(b,n))
                for representation,(start,end) in intervals.items():
                    assert 0<=start<end<=len(query[column])
                    identifiers=sorted((k for k in alignment if k!=sid),key=lambda k:(alignment[k],k))
                    seed=int(hashlib.sha256((family+'|'+role+'|'+query[column]).encode()).hexdigest()[:16],16)
                    random.Random(seed).shuffle(identifiers)
                    sequences=[cropped_a3m(alignment[sid],alignment[k],start,end) for k in [sid]+identifiers]
                    sequences=[sequences[0]]+[s for s in sequences[1:] if any(c.isupper() and c!='-' for c in s)]
                    assert sequences[0]==query[column][start:end]
                    assert all(len(''.join(c for c in s if not c.islower()))==end-start for s in sequences)
                    key=(family,role,query[column][start:end])
                    contents=canonical_msa.setdefault(key,msa_csv(sequences))
                    path=OUT/'msas'/(sid+'_'+role+'_'+representation+'.csv');path.write_text(contents)
                    msa_lookup[sid,role,representation]=path
    distance_path=ROOT/'results/art_expansion_20261105/step5/matching_pair_comparisons.tsv'
    with distance_path.open(newline='') as f:comparisons=list(csv.DictReader(f,delimiter='\t'))
    distance={frozenset([r['locus_a'],r['locus_b']]):1-(float(r['RT_core_identity'])+float(r['partner_identity']))/2 for r in comparisons}
    chosen=[min(typeI,key=lambda r:(int(r['RT_length_aa'])+int(r['partner_length_aa']),r['sample_id']))]
    while len(chosen)<config['domain_panel_size']:
        unused=[r for r in typeI if r not in chosen]
        def priority(r):
            novelty=min(distance.get(frozenset([r['locus_id'],s['locus_id']]),0) for s in chosen)
            new_cluster=all(r['partner_cluster70']!=s['partner_cluster70'] for s in chosen)
            new_source=all(r['source_group']!=s['source_group'] for s in chosen)
            return (new_source,new_cluster,novelty,-int(r['RT_length_aa'])-int(r['partner_length_aa']))
        chosen.append(max(unused,key=priority))
    table(OUT/'panels/domain_panel.csv',chosen)
    cases=[];request_index={}
    def add(stage,rt,partner,representation,condition,family,relation):
        entities=[];meta=[]
        if rt:
            n=len(rt['RT_protein_sequence']);b=int(rt['NTD_length_aa'])
            intervals=dict(full=(0,n),NTD_baseline=(0,b),NTD_minus20=(0,b-20),NTD_plus20=(0,b+20),core_only=(b,n))
            start,end=intervals[representation];sequence=rt['RT_protein_sequence'][start:end]
            entities.append(polymer('A',sequence,condition,msa_lookup.get((rt['sample_id'],'RT',representation))))
            meta.append(dict(chain_id='A',role='RT',sample_id=rt['sample_id'],start0=start,end0=end,source_NTD_end0=b))
        if partner:
            chain='B' if rt else 'A';sequence=partner['partner_protein_sequence']
            entities.append(polymer(chain,sequence,condition,msa_lookup.get((partner['sample_id'],'partner','full'))))
            meta.append(dict(chain_id=chain,role='partner',sample_id=partner['sample_id'],start0=0,end0=len(sequence)))
        payload=dict(polymers=entities,ligands=[],constraints=[],recycling_steps=config['recycling_steps'],sampling_steps=config['sampling_steps'],
            diffusion_samples=config['diffusion_samples'],step_scale=config['step_scale'],without_potentials=True,
            output_format='mmcif',concatenate_msas=False)
        digest=payload_hash(payload)
        request_path=OUT/'requests'/(digest+'.json')
        if not request_path.exists():write(request_path,payload)
        case_id=f'C{len(cases)+1:05d}'
        row=dict(case_id=case_id,stage=stage,family=family,representation=representation,condition=condition,
            RT_sample_id=rt['sample_id'] if rt else '',partner_sample_id=partner['sample_id'] if partner else '',
            source_pairing=relation,native_pair=bool(rt and partner and rt['sample_id']==partner['sample_id']),
            request_hash=digest,request_path=str(request_path.relative_to(OUT)),chain_metadata_json=json.dumps(meta),
            diffusion_samples=config['diffusion_samples'],RT_NTD_boundary_scope='operational_first_RVT1_envelope',confirmatory=False)
        cases.append(row)
        request_index.setdefault(digest,dict(request_hash=digest,request_path=row['request_path'],case_ids=[],stages=[]))
        request_index[digest]['case_ids'].append(case_id)
        if stage not in request_index[digest]['stages']:request_index[digest]['stages'].append(stage)
    # Native pairs go first, permitting early execution and calibration diagnostics
    # while all membership and alternatives remain fixed before scores are opened.
    for family,group in [('Type_I',typeI),('Type_II',typeII)]:
        for condition in config['conditions']:
            for rt in group:add('full_length_'+family,rt,rt,'full',condition,family,'native_source_pair')
            for rt in group:
                for partner in group:
                    if rt['sample_id']!=partner['sample_id']:add('full_length_'+family,rt,partner,'full',condition,family,'alternative_source_pair_unmeasured')
    for representation in config['domain_constructs']:
        for condition in config['domain_matrix_conditions']:
            for rt in chosen:
                for partner in chosen:add('domain_'+representation,rt,partner,representation,condition,'Type_I',
                    'native_source_pair' if rt['sample_id']==partner['sample_id'] else 'alternative_source_pair_unmeasured')
    for rt in chosen:
        for representation in ['full']+config['domain_constructs']:
            add('monomer',rt,None,representation,config['monomer_condition'],'Type_I','monomer')
        add('monomer',None,rt,'partner_only',config['monomer_condition'],'Type_I','monomer')
    n=config['cross_type_panel_size_each'];crossI=sorted(chosen,key=lambda r:r['sample_id'])[:n];crossII=typeII[:n]
    for condition in config['conditions']:
        for rt,partner in itertools.chain(itertools.product(crossI,crossII),itertools.product(crossII,crossI)):
            add('cross_type',rt,partner,'full',condition,'cross_type','cross_type_unmeasured')
    for r in review:
        for condition in config['review_conditions']:add('family_review_native',r,r,'full',condition,'Type_I_review','native_source_pair_family_unresolved')
    table(OUT/'experiments.csv',cases)
    requests=list(request_index.values())
    # Record request associations separately; content-identical cases are not
    # repredicted or presented as independent model samples.
    write(OUT/'request_index.json',requests)
    for name in ['art_bindings_common.py','art_bindings_api.py','art_bindings_prepare.py','art_bindings_run.py','art_bindings_report.py','art_bindings_tests.py','art_binding_api_credentials.py']:
        path=ROOT/'scripts'/name
        if path.is_file():shutil.copy2(path,OUT/'provenance/code'/name)
    shutil.copy2(source_config,OUT/'provenance/design_config.json')
    for path in [ROOT/'data/ART binding/api_pilot_20261105/report.json',ROOT/'data/ART binding/pilot_20261105/pilot_summary.json']:
        if path.exists():shutil.copy2(path,OUT/'provenance/prior_pilots'/(path.parent.name+'_'+path.name))
    artifacts=[source_config,source,distance_path,ROOT/'data/art_samples_2026-11-05/SHA256SUMS.json']
    protocol=dict(config=config,source_config_sha256=sha(source_config),prepared_utc=utc(),slurm_job_id=os.environ['SLURM_JOB_ID'],
        samples=dict(total=46,typeI_strict=18,typeII_complete=10,typeI_review=13,incomplete_excluded=5),
        cases=len(cases),unique_requests=len(requests),planned_models=len(requests)*config['diffusion_samples'],
        case_counts_by_stage=dict(Counter(c['stage'] for c in cases)),domain_panel=[r['sample_id'] for r in chosen],
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in artifacts},
        frozen_files={str(p.relative_to(OUT)):sha(p) for directory in ['inputs','panels','alignments','msas','requests','provenance/code'] for p in (OUT/directory).rglob('*') if p.is_file()},
        empirical_binding_calibration_complete=False,eligible_confirmatory_blocks=0,confirmatory=False)
    write(OUT/'protocol.json',protocol)
    write(OUT/'checkpoints/01_inputs_frozen.json',dict(status='complete',cases=len(cases),unique_requests=len(requests),planned_models=protocol['planned_models'],slurm_job_id=os.environ['SLURM_JOB_ID']))
    write(OUT/'status.json',dict(status='inputs_frozen',completed_requests=0,planned_requests=len(requests),planned_models=protocol['planned_models'],confirmatory=False))
    print(json.dumps({k:protocol[k] for k in ['samples','cases','unique_requests','planned_models','case_counts_by_stage','domain_panel']},indent=2),flush=True)


if __name__=='__main__':main()
