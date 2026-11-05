"""Execute frozen exploratory phases with bounded concurrency and checkpoints."""
import argparse
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
import json
import os
from pathlib import Path
import shutil
import signal
import threading
import time

from art_binding_api_credentials import credential
from art_bindings_common import ROOT,OUT,utc,sha,write,read_table
from art_bindings_api import run_request,archive_gateway_failures


def inventory(records,cached=None):
    statuses=Counter();models=0;stages=defaultdict(Counter)
    for record in records:
        path=OUT/'runs'/record['request_hash']/'state.json'
        state=cached.get(record['request_hash'],{}) if cached is not None else json.loads(path.read_text()) if path.exists() else {}
        status=state.get('status','not_started');statuses[status]+=1
        if status=='complete':models+=state.get('completed_models',0)
        for stage in record['stages']:stages[stage][status]+=1
    return dict(request_statuses=dict(statuses),completed_models=models,stages={k:dict(v) for k,v in stages.items()})


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['smoke','full'],required=True)
    parser.add_argument('--retry-gateway-failures',action='store_true')
    args=parser.parse_args()
    import fcntl
    lock=(OUT/'launcher.lock').open('w')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit('Another results_bindings launcher is already active.')
    protocol=json.loads((OUT/'protocol.json').read_text());config=protocol['config']
    for relative,digest in protocol['frozen_files'].items():assert sha(OUT/relative)==digest,relative
    for name,path_hash in protocol['source_sha256'].items():assert sha(ROOT/name)==path_hash,name
    key,provider,mode=credential(ROOT);assert provider=='nvidia_nim' and mode=='live'
    records=json.loads((OUT/'request_index.json').read_text());cases=read_table(OUT/'experiments.csv')
    if args.retry_gateway_failures:
        revisions=archive_gateway_failures(records)
        write(OUT/'checkpoints'/('transport_resume_'+os.environ['SLURM_JOB_ID']+'.json'),dict(
            status='recorded',scientific_payloads_unchanged=True,poll_seconds=300,retries=revisions,
            evidence='Gateway response without a result: either no returned ID or a terminal errored header; all attempts retained',updated_utc=utc()))
    cached={}
    for record in records:
        path=OUT/'runs'/record['request_hash']/'state.json'
        cached[record['request_hash']]=json.loads(path.read_text()) if path.exists() else dict(status='not_started')
    execution=OUT/'provenance/execution_code'/os.environ['SLURM_JOB_ID'];execution.mkdir(parents=True,exist_ok=True)
    code={}
    for name in ['art_bindings_common.py','art_bindings_api.py','art_bindings_run.py','art_binding_api_credentials.py']:
        source=ROOT/'scripts'/name;shutil.copy2(source,execution/name);code[name]=sha(source)
    write(execution/'context.json',dict(phase=args.phase,started_utc=utc(),code_sha256=code,
        cases_sha256=sha(OUT/'experiments.csv'),request_index_sha256=sha(OUT/'request_index.json')))
    if args.phase=='smoke':
        smoke=[]
        for family in ['Type_I','Type_II']:
            for condition in config['conditions']:
                smoke.append(next(c for c in cases if c['stage']=='full_length_'+family and c['condition']==condition and c['native_pair']=='True'))
        smoke.append(next(c for c in cases if c['stage']=='domain_NTD_baseline' and c['condition'].startswith('U_')))
        smoke.append(next(c for c in cases if c['stage']=='domain_core_only' and c['condition'].startswith('U_')))
        smoke.append(next(c for c in cases if c['stage']=='monomer' and c['RT_sample_id'] and c['representation']=='full'))
        smoke.append(next(c for c in cases if c['stage']=='monomer' and not c['RT_sample_id']))
        hashes={c['request_hash'] for c in smoke};selected=[r for r in records if r['request_hash'] in hashes]
    else:
        assert json.loads((OUT/'checkpoints/02_API_shapes_verified.json').read_text())['status']=='complete'
        selected=records
    errors=[];completed=0;consecutive_errors=0;start=time.monotonic();last_progress=0.0
    last_summary={}
    def progress(status):
        nonlocal last_progress,last_summary
        last_summary=dict(run_id=config['run_id'],status=status,phase=args.phase,updated_utc=utc(),
            slurm_job_id=os.environ['SLURM_JOB_ID'],planned_requests=len(records),planned_models=protocol['planned_models'],
            phase_selected_requests=len(selected),phase_returned_requests=completed,errors=errors,confirmatory=False,**inventory(records,cached))
        write(OUT/'status.json',last_summary)
        if last_summary.get('completed_models',0):write(OUT/'checkpoints/partial_progress.json',last_summary)
        for (stage,condition),members in case_groups.items():
            path=OUT/'checkpoints/stages'/(stage+'__'+condition+'.json')
            if not path.exists() and all(cached.get(h,{}).get('status')=='complete' for h in members):
                write(path,dict(status='complete',stage=stage,condition=condition,requests=len(members),
                    model_samples=sum(cached[h].get('completed_models',0) for h in members),completed_utc=utc(),
                    slurm_job_id=os.environ['SLURM_JOB_ID'],confirmatory=False))
        last_progress=time.monotonic()
    case_groups=defaultdict(set)
    for case in cases:case_groups[case['stage'],case['condition']].add(case['request_hash'])
    progress('running')
    iterator=iter(selected);pending={};stop=False
    def request_stop(signum,frame):
        nonlocal stop
        stop=True
        print('Graceful stop requested; finishing in-flight invocations.',flush=True)
    signal.signal(signal.SIGINT,request_stop);signal.signal(signal.SIGTERM,request_stop)
    with ThreadPoolExecutor(max_workers=config['parallel_requests']) as executor:
        def fill():
            while not stop and len(pending)<config['parallel_requests']:
                try:record=next(iterator)
                except StopIteration:break
                if cached.get(record['request_hash'],{}).get('status')!='complete':cached[record['request_hash']]=dict(status='running')
                pending[executor.submit(run_request,record,key,config)]=record
        fill()
        while pending:
            ready,_=wait(pending,timeout=30,return_when=FIRST_COMPLETED)
            for future in ready:
                record=pending.pop(future)
                try:
                    result=future.result();completed+=1;consecutive_errors=0;cached[record['request_hash']]=result
                    print(f'COMPLETE {completed}/{len(selected)} {record["request_hash"][:12]} {record["stages"][0]}',flush=True)
                except Exception as error:
                    message=str(error).replace(key,'[REDACTED]');errors.append(dict(request_hash=record['request_hash'],error=message));consecutive_errors+=1
                    path=OUT/'runs'/record['request_hash']/'state.json'
                    cached[record['request_hash']]=json.loads(path.read_text()) if path.exists() else dict(status='failed',error=message)
                    print('FAILED '+record['request_hash'][:12]+' '+message,flush=True)
                    if args.phase=='smoke' or consecutive_errors>=12 or any('HTTP '+str(c)+' ' in message for c in [401,402,403]):stop=True
            if ready or time.monotonic()-last_progress>=60:progress('running')
            fill()
    result_inventory=inventory(records,cached)
    successful=not errors and completed==len(selected)
    progress('complete' if successful and args.phase=='full' else 'smoke_complete' if successful else 'incomplete')
    checkpoint='02_API_shapes_verified.json' if args.phase=='smoke' else '03_predictions_complete.json'
    write(OUT/'checkpoints'/checkpoint,dict(status='complete' if successful else 'incomplete',phase=args.phase,
        selected_requests=len(selected),returned_requests=completed,errors=errors,seconds=time.monotonic()-start,
        slurm_job_id=os.environ['SLURM_JOB_ID'],**result_inventory))
    print(json.dumps({k:last_summary[k] for k in ['status','phase','planned_requests','planned_models','request_statuses','completed_models']},indent=2),flush=True)
    if not successful:raise SystemExit('Phase incomplete; failure states and pending request IDs preserved.')


if __name__=='__main__':main()
