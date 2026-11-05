"""Export all completed exploratory models, balanced summaries and figures."""
import argparse
from collections import Counter,defaultdict
import csv
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import statistics

import numpy as np

from art_bindings_common import ROOT,OUT,utc,sha,write,table,read_table


def numeric(value):return isinstance(value,(int,float)) and math.isfinite(value)


def median(values):
    values=[v for v in values if numeric(v)]
    return float(statistics.median(values)) if values else None


def interface_row(case,model,folder):
    row={k:case[k] for k in ['case_id','stage','family','representation','condition','RT_sample_id','partner_sample_id','source_pairing','native_pair','request_hash']}
    scores=model['scores']
    row.update(model_index=model['model_index'],protein_iptm=scores['protein_iptm_scores'],iptm=scores['iptm_scores'],
        ptm=scores['ptm_scores'],confidence=scores['confidence_scores'],complex_plddt=scores['complex_plddt_scores'],
        structure_file=str((folder/model['structure_file']).relative_to(OUT)),structure_sha256=model['structure_sha256'],
        contacts_4A=model['contacts_by_cutoff']['4'],contacts_5A=model['contacts_by_cutoff']['5'],contacts_8A=model['contacts_by_cutoff']['8'],
        chain0_mean_CA_plddt=model['chain_summary'][0]['mean_CA_plddt'],
        chain1_mean_CA_plddt=model['chain_summary'][1]['mean_CA_plddt'] if len(model['chain_summary'])==2 else None,
        fixed_seed_available=False,confirmatory=False)
    meta=json.loads(case['chain_metadata_json'])
    arrays=np.load(folder/model['coordinates_file'])
    if meta[0]['role']=='RT':
        r=meta[0];n=len(arrays['chain0_CA']);local_boundary=max(0,min(n,r['source_NTD_end0']-r['start0']))
        confidence=arrays['chain0_plddt']
        row.update(NTD_residues_in_construct=local_boundary,core_residues_in_construct=n-local_boundary,
            NTD_mean_CA_plddt=float(confidence[:local_boundary].mean()) if local_boundary else None,
            core_mean_CA_plddt=float(confidence[local_boundary:].mean()) if n>local_boundary else None)
        if model['contacts_file']:
            contacts=read_table(folder/model['contacts_file']);ntd_contacts=[];core_contacts=[];confident=[]
            for c in contacts:
                i=int(c['chain0_residue0']);j=int(c['chain1_residue0']);d=float(c['min_heavy_atom_distance_A'])
                if d>5:continue
                (ntd_contacts if i<local_boundary else core_contacts).append((i,j))
                if confidence[i]>=70 and arrays['chain1_plddt'][j]>=70:confident.append((i,j))
            all_contacts=len(ntd_contacts)+len(core_contacts)
            ntd_unique=len({i for i,j in ntd_contacts});core_unique=len({i for i,j in core_contacts})
            row.update(NTD_contact_pairs_5A=len(ntd_contacts),core_contact_pairs_5A=len(core_contacts),
                NTD_fraction_of_contact_pairs=len(ntd_contacts)/all_contacts if all_contacts else None,
                NTD_fraction_residues_contacting=ntd_unique/local_boundary if local_boundary else None,
                core_fraction_residues_contacting=core_unique/(n-local_boundary) if n>local_boundary else None,
                confident_contact_pairs_5A=len(confident))
    return row


def rmsd(reference,moving):
    first=reference-reference.mean(axis=0);second=moving-moving.mean(axis=0)
    u,s,v=np.linalg.svd(second.T@first);rotation=u@v
    if np.linalg.det(rotation)<0:u[:,-1]*=-1;rotation=u@v
    return float(np.sqrt(np.mean(np.sum((second@rotation-first)**2,axis=1))))


def monomer_integrity(cases,normalized):
    groups=defaultdict(dict)
    for c in cases:
        if c['stage']=='monomer' and c['RT_sample_id'] and c['request_hash'] in normalized:
            groups[c['RT_sample_id']][c['representation']]=c
    comparisons=[]
    for sid,representations in groups.items():
        if 'full' not in representations:continue
        full=representations['full'];full_folder=OUT/'runs'/full['request_hash']
        for representation,case in representations.items():
            if representation=='full':continue
            meta=json.loads(case['chain_metadata_json'])[0];a,b=meta['start0'],meta['end0'];folder=OUT/'runs'/case['request_hash']
            for m1,m2 in itertools.product(normalized[full['request_hash']]['models'],normalized[case['request_hash']]['models']):
                r=np.load(full_folder/m1['coordinates_file']);q=np.load(folder/m2['coordinates_file'])
                reference=r['chain0_CA'][a:b];moving=q['chain0_CA'];assert reference.shape==moving.shape
                selected=(r['chain0_plddt'][a:b]>=70)&(q['chain0_plddt']>=70)
                comparisons.append(dict(sample_id=sid,representation=representation,full_model=m1['model_index'],construct_model=m2['model_index'],
                    compared_residues=len(reference),CA_RMSD_A=rmsd(reference,moving),confident_matched_residues=int(selected.sum()),
                    confident_CA_RMSD_A=rmsd(reference[selected],moving[selected]) if selected.sum()>=20 else None,
                    full_region_mean_plddt=float(r['chain0_plddt'][a:b].mean()),construct_mean_plddt=float(q['chain0_plddt'].mean()),
                    boundary_scope='operational_RVT1_envelope',experimental_reference=False))
    return comparisons


def plot_matrix(plt,matrix,labels,path,title):
    fig,axis=plt.subplots(figsize=(max(6,len(labels)*.43),max(5,len(labels)*.40)))
    cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#dddddd')
    display=axis.imshow(np.ma.masked_invalid(matrix),vmin=0,vmax=1,cmap=cmap)
    axis.set_xticks(range(len(labels)),labels=labels,rotation=90,fontsize=7)
    axis.set_yticks(range(len(labels)),labels=labels,fontsize=7)
    axis.set_xlabel('Partner source');axis.set_ylabel('RT source');axis.set_title(title,fontsize=10)
    from matplotlib.patches import Rectangle
    for i in range(len(labels)):axis.add_patch(Rectangle((i-.5,i-.5),1,1,fill=False,edgecolor='white',linewidth=.8))
    fig.colorbar(display,ax=axis,label='Median protein ipTM (model confidence)')
    fig.text(.02,.01,'Outlined diagonal: native source pairs. Missing results are grey. Exploratory predictions.',fontsize=7)
    fig.tight_layout(rect=[0,.03,1,1])
    for ext in ['png','pdf','svg']:fig.savefig(path.with_suffix('.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    parser=argparse.ArgumentParser();parser.add_argument('--final',action='store_true');args=parser.parse_args()
    protocol=json.loads((OUT/'protocol.json').read_text());cases=read_table(OUT/'experiments.csv');index=json.loads((OUT/'request_index.json').read_text())
    normalized={};states={};failures=[]
    for request in index:
        folder=OUT/'runs'/request['request_hash'];path=folder/'state.json'
        state=json.loads(path.read_text()) if path.exists() else dict(status='not_started')
        states[request['request_hash']]=state
        if state['status']=='complete' and (folder/'normalized.json').exists():
            data=json.loads((folder/'normalized.json').read_text());assert data['verified']
            for model in data['models']:assert sha(folder/model['structure_file'])==model['structure_sha256']
            normalized[request['request_hash']]=data
        elif state['status'] in ['failed','submission_uncertain']:
            failures.append(dict(request_hash=request['request_hash'],status=state['status'],error=state.get('error',''),request_id=state.get('request_id')))
    model_rows=[];pair_rows=[]
    for case in cases:
        data=normalized.get(case['request_hash'])
        if not data:continue
        group=[interface_row(case,m,OUT/'runs'/case['request_hash']) for m in data['models']];model_rows.extend(group)
        row={k:case[k] for k in ['case_id','stage','family','representation','condition','RT_sample_id','partner_sample_id','source_pairing','native_pair','request_hash']}
        row.update(model_samples=len(group),distinct_structure_hashes=len({m['structure_sha256'] for m in group}),
            median_protein_iptm=median(m['protein_iptm'] for m in group),min_protein_iptm=min((m['protein_iptm'] for m in group if numeric(m['protein_iptm'])),default=None),
            max_protein_iptm=max((m['protein_iptm'] for m in group if numeric(m['protein_iptm'])),default=None),
            median_ptm=median(m['ptm'] for m in group),median_confidence=median(m['confidence'] for m in group),
            median_contacts_5A=median(m['contacts_5A'] for m in group),
            median_NTD_fraction_contacts=median(m.get('NTD_fraction_of_contact_pairs') for m in group),
            median_confident_contacts_5A=median(m.get('confident_contact_pairs_5A') for m in group),confirmatory=False)
        pair_rows.append(row)
    # Dict field union preserves optional domain/monomer values without turning
    # absent values into zero or dropping incomplete analysis groups.
    fields=list(dict.fromkeys(k for r in model_rows for k in r));table(OUT/'results/model_scores.csv',model_rows,fields or None)
    table(OUT/'results/pair_medians.csv',pair_rows);table(OUT/'results/failures.csv',failures)
    attempt_rows=[]
    for request in index:
        folder=OUT/'runs'/request['request_hash']
        for path in sorted((folder/'attempts').glob('*/retry_decision.json')):
            attempt=json.loads(path.read_text())
            attempt_rows.append(dict(attempt,current_status=states[request['request_hash']]['status']))
    table(OUT/'results/retried_gateway_attempts.csv',attempt_rows)
    table(OUT/'results/monomer_domain_integrity.csv',monomer_integrity(cases,normalized))
    contrast_rows=[];row_contrasts=[];matrices={}
    domain_ids=[r['sample_id'] for r in read_table(OUT/'panels/domain_panel.csv')]
    for (stage,condition),rows in itertools.groupby(sorted(pair_rows,key=lambda r:(r['stage'],r['condition'])),key=lambda r:(r['stage'],r['condition'])):
        rows=list(rows)
        if stage not in ['full_length_Type_I','full_length_Type_II'] and not stage.startswith('domain_'):continue
        if stage=='full_length_Type_I':labels=[r['sample_id'] for r in read_table(OUT/'panels/typeI_strict.csv')]
        elif stage=='full_length_Type_II':labels=[r['sample_id'] for r in read_table(OUT/'panels/typeII_complete.csv')]
        else:labels=domain_ids
        matrix=np.full((len(labels),len(labels)),np.nan);lookup={label:i for i,label in enumerate(labels)}
        for r in rows:
            if numeric(r['median_protein_iptm']):matrix[lookup[r['RT_sample_id']],lookup[r['partner_sample_id']]]=r['median_protein_iptm']
        matrices[stage,condition]=(labels,matrix)
        table(OUT/'results'/(stage+'__'+condition+'__matrix.csv'),[dict(RT_sample_id=label,**{p:None if np.isnan(matrix[i,j]) else matrix[i,j] for j,p in enumerate(labels)}) for i,label in enumerate(labels)])
        complete=bool(np.isfinite(matrix).all())
        native=np.diag(matrix);alternatives=matrix[~np.eye(len(labels),dtype=bool)]
        contrast_rows.append(dict(stage=stage,condition=condition,loci=len(labels),completed_cells=int(np.isfinite(matrix).sum()),planned_cells=matrix.size,
            complete_matrix=complete,native_mean=float(native.mean()) if complete else None,alternative_mean=float(alternatives.mean()) if complete else None,
            descriptive_diagonal_minus_offdiagonal=float(native.mean()-alternatives.mean()) if complete else None,
            independent_blocks_established=0,pvalue=None,confirmatory=False))
        for i,label in enumerate(labels):
            others=np.delete(matrix[i],i)
            if np.isfinite(matrix[i,i]) and np.isfinite(others).all():
                row_contrasts.append(dict(stage=stage,condition=condition,RT_sample_id=label,
                    native_score=float(matrix[i,i]),alternative_mean=float(others.mean()),row_contrast=float(matrix[i,i]-others.mean())))
    # Domain representations are compared with the same six-locus full-length
    # submatrix; never compare their six-locus contrast with an eighteen-locus one.
    matched=[]
    for condition in protocol['config']['conditions']:
        base=matrices.get(('full_length_Type_I',condition))
        if not base:continue
        labels,matrix=base;indices=[labels.index(s) for s in domain_ids];matrix=matrix[np.ix_(indices,indices)]
        for representation,selection in [('full',matrix)]+[(name,matrices.get(('domain_'+name,condition),(None,None))[1]) for name in protocol['config']['domain_constructs']]:
            if selection is None or not np.isfinite(selection).all():continue
            diag=np.diag(selection);off=selection[~np.eye(len(domain_ids),dtype=bool)]
            matched.append(dict(condition=condition,representation=representation,loci=len(domain_ids),
                native_mean=float(diag.mean()),alternative_mean=float(off.mean()),contrast=float(diag.mean()-off.mean()),confirmatory=False))
    table(OUT/'results/matrix_contrasts.csv',contrast_rows);table(OUT/'results/RT_row_contrasts.csv',row_contrasts)
    table(OUT/'results/matched_domain_panel_contrasts.csv',matched)
    plots=False
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for (stage,condition),(labels,matrix) in matrices.items():
            plot_matrix(plt,matrix,labels,OUT/'figures'/(stage+'__'+condition),stage+' · '+('query-only' if condition.startswith('S_') else 'requested unpaired MSA'))
        if matched:
            fig,axis=plt.subplots(figsize=(8,4));names=['full']+protocol['config']['domain_constructs']
            for offset,condition in enumerate(protocol['config']['conditions']):
                values={r['representation']:r['contrast'] for r in matched if r['condition']==condition}
                axis.plot(range(len(names)),[values.get(name,np.nan) for name in names],marker='o',label=condition)
            axis.axhline(0,color='grey',linewidth=.8);axis.set_xticks(range(len(names)),names,rotation=20,ha='right')
            axis.set_ylabel('Native minus alternative median-score means');axis.set_title('Matched six-locus domain panel: exploratory contrast');axis.legend(fontsize=7)
            fig.tight_layout()
            for extension in ['png','pdf','svg']:fig.savefig(OUT/'figures'/('matched_domain_contrasts.'+extension),dpi=300,bbox_inches='tight')
            plt.close(fig)
        plots=True
    except ImportError:pass
    completed=len(normalized);models=sum(len(v['models']) for v in normalized.values())
    status_counts=dict(Counter(s['status'] for s in states.values()))
    report=dict(run_id=protocol['config']['run_id'],updated_utc=utc(),completed_requests=completed,planned_requests=len(index),
        completed_model_samples=models,planned_model_samples=protocol['planned_models'],completed_case_rows=len(pair_rows),request_statuses=status_counts,
        complete=completed==len(index),figures_generated=plots,biological_binding_validated=False,empirical_specificity_calibration_complete=False,
        confirmed_independent_blocks=0,confirmatory=False,matrix_summaries=contrast_rows,
        limitations=protocol['config']['limitations'],slurm_job_id=os.environ['SLURM_JOB_ID'])
    write(OUT/'results/summary.json',report)
    planned_groups=Counter((c['stage'],c['condition']) for c in cases)
    finished_groups=Counter((c['stage'],c['condition']) for c in pair_rows)
    lines=['# ART binding prediction results', '',
        '**Status: '+('complete prediction batch' if report['complete'] else 'partial prediction batch')+'** · '+report['updated_utc'], '',
        f"{completed:,} of {len(index):,} API requests returned {models:,} of {protocol['planned_models']:,} planned structure models. "
        'Every counted model passed input-chain sequence verification and an immutable structure checksum check.', '',
        'These are predicted structures and model confidence scores. Physical binding and native-pair specificity have not been established.', '',
        '| Experiment | Input condition | Completed / planned requests |', '|---|---|---:|']
    for (stage,condition),planned in planned_groups.items():
        lines.append(f'| {stage} | {condition} | {finished_groups[stage,condition]} / {planned} |')
    lines.extend(['', '## Available outputs', '',
        '- [Scores and contacts for each model](model_scores.csv)',
        '- [Three-model medians for each pairing](pair_medians.csv)',
        '- [Matrix coverage and descriptive contrasts](matrix_contrasts.csv)',
        '- [Matched six-locus domain comparisons](matched_domain_panel_contrasts.csv)',
        '- [Monomer domain integrity checks](monomer_domain_integrity.csv)',
        '- [Current failures](failures.csv) and [archived gateway retry attempts](retried_gateway_attempts.csv)',
        '- [Figures in PNG, PDF and SVG](../figures/)', '',
        '## Interpretation', '',
        'A native-versus-alternative contrast is exported only after its complete comparison matrix has returned. '
        'Missing cells remain missing. Low interface confidence or a null contrast cannot establish absence of biological binding. '
        'Three stochastic models per cell describe prediction variability; they are not three biological samples.', '',
        'The U-requested arm supplies independent per-chain homolog alignments, but the hosted service does not expose '
        'its processed features. It remains an exploratory MSA sensitivity analysis. The NTD cuts use operational RVT1 boundaries.', '',
        'Zero eligible independent three-locus blocks are established. Empirical specificity calibration, verified NART controls, '
        'and the RNA/ternary extension remain unresolved; no confirmatory inference is reported.', ''])
    (OUT/'results/REPORT.md').write_text('\n'.join(lines))
    table(OUT/'results/deferred_arms.csv',[
        dict(arm='verified_strict_U_condition',reason='Hosted processed features unavailable; U-requested is exploratory'),
        dict(arm='NART_controls',reason='No verified accession-linked NART controls available'),
        dict(arm='empirical_specificity_calibration',reason='Two experimentally characterized comparison families have not been assembled and validated'),
        dict(arm='RNA_and_ternary_extension',reason='Optional arm; RNA boundaries, reference-complex calibration and expanded array validation remain unresolved'),
        dict(arm='confirmatory_inference',reason='Zero eligible matched evolutionary blocks; calibration and power gates unmet')])
    if args.final:
        source=Path(__file__);shutil.copy2(source,OUT/'provenance'/'report_code.py')
        hashes={str(p.relative_to(OUT)):sha(p) for directory in ['inputs','panels','alignments','msas','requests','results','figures','provenance'] for p in (OUT/directory).rglob('*') if p.is_file()}
        write(OUT/'report_artifact_sha256.json',hashes)
        write(OUT/'checkpoints/04_analysis_exported.json',dict(status='complete' if report['complete'] else 'partial',completed_requests=completed,
            completed_models=models,figures_generated=plots,artifacts_hashed=len(hashes),slurm_job_id=os.environ['SLURM_JOB_ID']))
    print(json.dumps({k:report[k] for k in ['completed_requests','planned_requests','completed_model_samples','planned_model_samples','complete','figures_generated','request_statuses']},indent=2),flush=True)


if __name__=='__main__':main()
