"""Execute the initial bounded development pilot and retain cell-level state."""
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/ART binding/pilot_20261105'


def write(path,value):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    manifest=json.loads((OUT/'manifest.json').read_text())
    provenance=json.loads((ROOT/'data/ART binding/env/boltz_provenance.json').read_text());assert provenance['passed']
    for relative,expected in manifest['hashes'].items():
        assert hashlib.sha256((OUT/relative).read_bytes()).hexdigest()==expected,relative
    for name in ['art_binding_boltz_entry.py','art_binding_msa_control.py']:
        assert hashlib.sha256((ROOT/'scripts'/name).read_bytes()).hexdigest()==hashlib.sha256((OUT/'code'/name).read_bytes()).hexdigest(),name
    config=manifest['config'];states=[];started=time.monotonic()
    for cell in manifest['cells']:
        if time.monotonic()-started>9000:break
        folder=OUT/'runs'/(cell['cell_id']+'_seed'+str(cell['seed']));folder.mkdir(parents=True,exist_ok=True)
        status_path=folder/'status.json'
        if status_path.exists():
            previous=json.loads(status_path.read_text())
            if previous.get('status')=='complete':states.append(previous);continue
        path=ROOT/cell['input_path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==cell['input_sha256']
        args=[sys.executable,str(ROOT/'scripts/art_binding_boltz_entry.py'),'predict',str(path),
            '--out_dir',str(folder),'--cache',str(ROOT/'.tools/boltz_cache'),'--model','boltz2',
            '--checkpoint',str(ROOT/'.tools/boltz_cache/boltz2_conf.ckpt'),
            '--accelerator','gpu','--devices','1','--seed',str(cell['seed']),
            '--recycling_steps',str(config['recycling_steps']),'--sampling_steps',str(config['sampling_steps']),
            '--diffusion_samples','1','--max_parallel_samples','1','--max_msa_seqs','64',
            '--num_workers','0','--preprocessing-threads','2','--no_kernels','--write_full_pae']
        env=dict(os.environ,ART_BINDING_CONDITION=cell['condition'],ART_BINDING_AUDIT_DIR=str(folder/'audit'),
            WANDB_MODE='disabled',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1')
        state=dict(cell,slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),status='running',
            started_utc=datetime.now(timezone.utc).isoformat(),command=args,attempts=[])
        write(status_path,state);print('RUNNING '+cell['cell_id'],flush=True)
        for attempt in range(2):
            start=time.monotonic()
            try:
                with (folder/f'attempt_{attempt+1}.log').open('w') as log:
                    process=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=1800)
                code=process.returncode;error=None
            except subprocess.TimeoutExpired:code=-1;error='1800-second cell timeout'
            state['attempts'].append(dict(returncode=code,error=error,seconds=time.monotonic()-start))
            confidence=list(folder.glob('boltz_results_*/predictions/*/confidence_*_model_0.json'))
            audits=list((folder/'audit').glob('feature_audit_*.json'))
            if code==0 and len(confidence)==1 and audits:
                state.update(status='complete',confidence_file=str(confidence[0].relative_to(ROOT)))
                break
            state['status']='failed';write(status_path,state)
            # Missing or rejected input features need correction before further
            # GPU work; a technical retry cannot bypass their safeguards.
            if not audits:break
        state['ended_utc']=datetime.now(timezone.utc).isoformat();states.append(state);write(status_path,state)
        write(OUT/'run_status.json',dict(status='running',cells=states,planned_cells=len(manifest['cells']),confirmatory=False))
        if state['status']!='complete':raise RuntimeError('Pilot cell failed; preserved evidence at '+str(folder))
    complete=sum(s['status']=='complete' for s in states)
    result=dict(status='complete' if complete==len(manifest['cells']) else 'partial_time_budget',
        cells=states,planned_cells=len(manifest['cells']),completed_cells=complete,confirmatory=False,
        API_key_required=False,slurm_job_id=os.environ['SLURM_JOB_ID'])
    write(OUT/'run_status.json',result);print(json.dumps(dict(status=result['status'],completed_cells=complete),indent=2),flush=True)


if __name__=='__main__':main()
