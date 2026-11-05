"""Archive the failed wire-format revision and freeze a semantic-equivalent fix."""
from collections import Counter
import csv
import io
import json
import os
from pathlib import Path
import shutil

from art_bindings_common import ROOT,OUT,utc,sha,write,table,read_table,payload_hash


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    protocol=json.loads((OUT/'protocol.json').read_text())
    history=OUT/'provenance/history/schema_header_revision1';history.mkdir(parents=True,exist_ok=False)
    for path in [OUT/'protocol.json',OUT/'experiments.csv',OUT/'request_index.json',OUT/'status.json',OUT/'checkpoints/02_API_shapes_verified.json',OUT/'provenance/design_config.json']:
        if path.exists():shutil.copy2(path,history/(str(path.relative_to(OUT)).replace('/','__')))
    shutil.copytree(OUT/'provenance/code',history/'code')
    cases=read_table(OUT/'experiments.csv');new_index={};changed=0
    for case in cases:
        path=OUT/case['request_path'];payload=json.loads(path.read_text())
        for p in payload['polymers']:
            record=p['msa']['controlled_input'].get('csv')
            if record:
                rows=list(csv.DictReader(io.StringIO(record['alignment'])))
                assert all(not r['key'] for r in rows)
                stream=io.StringIO();writer=csv.writer(stream,lineterminator='\n');writer.writerow(['key','sequence'])
                for r in rows:writer.writerow(['',r['sequence']])
                record['alignment']=stream.getvalue()
        digest=payload_hash(payload)
        changed+=digest!=case['request_hash']
        case['request_hash']=digest;case['request_path']='requests/'+digest+'.json'
        if not (OUT/case['request_path']).exists():write(OUT/case['request_path'],payload)
        row=new_index.setdefault(digest,dict(request_hash=digest,request_path=case['request_path'],case_ids=[],stages=[]))
        row['case_ids'].append(case['case_id'])
        if case['stage'] not in row['stages']:row['stages'].append(case['stage'])
    table(OUT/'experiments.csv',cases);write(OUT/'request_index.json',list(new_index.values()))
    source_config=ROOT/'config/art_bindings_exploratory_20261105.json'
    cfg=json.loads(source_config.read_text());assert cfg['parallel_requests']==2
    shutil.copy2(source_config,OUT/'provenance/design_config.json')
    protocol['config']=cfg;protocol['source_config_sha256']=sha(source_config)
    protocol['source_sha256'][str(source_config.relative_to(ROOT))]=sha(source_config)
    amendment=dict(completed_utc=utc(),reason='NVIDIA hosted validation requires CSV key,sequence order; request values unchanged. Reduce concurrency after rate-limit responses.',
        changed_payloads=changed,sampling_unchanged=True,panel_membership_unchanged=True,original_snapshot=str(history.relative_to(OUT)),slurm_job_id=os.environ['SLURM_JOB_ID'])
    protocol.setdefault('amendments',[]).append(amendment)
    for name in ['art_bindings_common.py','art_bindings_api.py','art_bindings_prepare.py','art_bindings_run.py','art_bindings_report.py','art_bindings_tests.py','art_binding_api_credentials.py','art_bindings_schema_amend.py']:
        source=ROOT/'scripts'/name
        if source.exists():shutil.copy2(source,OUT/'provenance/code'/name)
    protocol['frozen_files']={str(p.relative_to(OUT)):sha(p) for directory in ['inputs','panels','alignments','msas','requests','provenance/code'] for p in (OUT/directory).rglob('*') if p.is_file()}
    protocol['frozen_files'].update({name:sha(OUT/name) for name in ['experiments.csv','request_index.json']})
    write(OUT/'protocol.json',protocol);write(OUT/'checkpoints/01b_request_schema_corrected.json',dict(status='complete',**amendment))
    print(json.dumps(amendment,indent=2),flush=True)


if __name__=='__main__':main()
