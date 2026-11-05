"""Measure a bounded CPU fallback without changing registered GPU experiments."""
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import socket
import subprocess
import sys
import time


def write(path,value):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    root=Path(__file__).resolve().parents[1]
    out=root/'data/ART binding/cpu_probe_20261105'/os.environ['SLURM_JOB_ID'];out.mkdir(parents=True,exist_ok=True)
    pilot=root/'data/ART binding/pilot_20261105'
    manifest=json.loads((pilot/'manifest.json').read_text())
    cell=next(c for c in manifest['cells'] if c['cell_id']=='S_r1_p1')
    source=root/cell['input_path'];assert hashlib.sha256(source.read_bytes()).hexdigest()==cell['input_sha256']
    path=out/'S_r1_p1.yaml';shutil.copy2(source,path)
    args=[sys.executable,str(root/'scripts/art_binding_cpu_entry.py'),'predict',str(path),
        '--out_dir',str(out),'--cache',str(root/'.tools/boltz_cache'),'--model','boltz2',
        '--checkpoint',str(root/'.tools/boltz_cache/boltz2_conf.ckpt'),'--accelerator','cpu','--devices','1',
        '--seed',str(cell['seed']),'--recycling_steps','1','--sampling_steps','10',
        '--diffusion_samples','1','--max_parallel_samples','1','--max_msa_seqs','64',
        '--num_workers','0','--preprocessing-threads','2','--no_kernels','--write_full_pae']
    env=dict(os.environ,ART_BINDING_CONDITION='S',ART_BINDING_AUDIT_DIR=str(out/'audit'),WANDB_MODE='disabled',PYTHONDONTWRITEBYTECODE='1')
    result=dict(status='running',purpose='CPU_timing_and_compatibility_only',
        started_utc=datetime.now(timezone.utc).isoformat(),slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
        allocated_CPUs=int(os.environ['SLURM_CPUS_PER_TASK']),accelerator='cpu',precision='32-true',
        sequence_lengths=[cell['RT_length_aa'],cell['partner_length_aa']],
        sampling_steps=10,recycling_steps=1,protocol_sampling_steps=200,protocol_recycling_steps=3,
        timeout_seconds=1200,command=args,input_sha256=cell['input_sha256'],API_key_required=False,
        comparable_to_registered_GPU_results=False,
        limitation='Reduced-settings timing probe in FP32; never combine its scores with the registered GPU matrix.')
    write(out/'status.json',result)
    for name in ['art_binding_cpu_entry.py','art_binding_msa_control.py','art_binding_cpu_probe.py']:
        shutil.copy2(root/'scripts'/name,out/name)
    start=time.monotonic()
    try:
        with (out/'prediction.log').open('w') as log:
            process=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=1200)
        result['returncode']=process.returncode
        confidence=list(out.glob('boltz_results_*/predictions/*/confidence_*_model_0.json'))
        result['status']='complete' if process.returncode==0 and len(confidence)==1 else 'failed'
        if confidence:result['confidence_file']=str(confidence[0].relative_to(root))
    except subprocess.TimeoutExpired:
        result['status']='timing_budget_exhausted';result['returncode']=None
    result.update(seconds=time.monotonic()-start,ended_utc=datetime.now(timezone.utc).isoformat(),
        peak_child_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    write(out/'status.json',result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
