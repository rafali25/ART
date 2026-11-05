"""Run the pinned model with audited query-only or explicit unpaired features."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import torch

from art_binding_msa_control import separate_chains,audit_msa


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    condition=os.environ['ART_BINDING_CONDITION'];out=Path(os.environ['ART_BINDING_AUDIT_DIR'])
    out.mkdir(parents=True,exist_ok=True)
    from boltz.data import const
    import boltz.data.feature.featurizerv2 as features
    import boltz.main as boltz_main
    original=features.construct_paired_msa;calls=0
    def controlled(*args,**kwargs):
        nonlocal calls
        data=kwargs.get('data',args[0] if args else None)
        msa,deletions,paired=original(*args,**kwargs)
        chains=np.asarray(data.tokens['asym_id']);gap=const.token_ids['-']
        raw=[x.detach().cpu().numpy() for x in [msa,deletions,paired]]
        if condition=='U':
            changed=separate_chains(*raw,chains,gap)
            msa,deletions,paired=[torch.as_tensor(a,dtype=t.dtype,device=t.device) for a,t in zip(changed,[msa,deletions,paired])]
        checked=[x.detach().cpu().numpy() for x in [msa,deletions,paired]]
        audit=audit_msa(checked[0],checked[2],chains,gap,condition)
        assert audit['rows']<=64,'MSA exceeds the registered row cap after separation'
        audit['native_rows_before_control']=int(raw[0].shape[1])
        audit['adapter_sha256']=hashlib.sha256(Path(__file__).with_name('art_binding_msa_control.py').read_bytes()).hexdigest()
        np.savez_compressed(out/f'msa_features_{calls}.npz',msa=checked[0],deletion=checked[1],paired=checked[2],chain_ids=chains,
            native_msa=raw[0],native_paired=raw[2])
        (out/f'feature_audit_{calls}.json').write_text(json.dumps(audit,indent=2)+'\n')
        calls+=1
        return msa,deletions,paired
    features.construct_paired_msa=controlled
    assert torch.cuda.is_available() and torch.cuda.device_count()==1,'Require exactly one Slurm-visible GPU'
    torch.cuda.reset_peak_memory_stats();start=time.monotonic();succeeded=False
    try:
        boltz_main.cli.main(args=sys.argv[1:],standalone_mode=False)
        assert calls>0,'No actual MSA features were audited'
        succeeded=True
    finally:
        metrics=dict(success=succeeded,seconds=time.monotonic()-start,feature_audits=calls,
            peak_GPU_allocated_bytes=torch.cuda.max_memory_allocated(),peak_GPU_reserved_bytes=torch.cuda.max_memory_reserved(),
            GPU=torch.cuda.get_device_name(0),slurm_job_id=os.environ['SLURM_JOB_ID'])
        (out/'execution_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')


if __name__=='__main__':main()
