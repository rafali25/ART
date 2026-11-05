"""Audited CPU-only timing probe; separate from the frozen GPU pilot."""
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
import torch

from art_binding_msa_control import separate_chains,audit_msa


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use a Slurm CPU job.')
    if '--accelerator' not in sys.argv or sys.argv[sys.argv.index('--accelerator')+1]!='cpu':
        raise SystemExit('This entry point only permits CPU execution.')
    threads=int(os.environ['SLURM_CPUS_PER_TASK']);torch.set_num_threads(threads);torch.set_num_interop_threads(2)
    condition=os.environ['ART_BINDING_CONDITION'];out=Path(os.environ['ART_BINDING_AUDIT_DIR']);out.mkdir(parents=True,exist_ok=True)
    from boltz.data import const
    import boltz.data.feature.featurizerv2 as features
    import boltz.main as boltz_main
    native_construct=features.construct_paired_msa;calls=0
    def controlled(*args,**kwargs):
        nonlocal calls
        data=kwargs.get('data',args[0] if args else None)
        msa,deletion,paired=native_construct(*args,**kwargs)
        chains=np.asarray(data.tokens['asym_id']);gap=const.token_ids['-']
        arrays=[t.detach().cpu().numpy() for t in [msa,deletion,paired]]
        if condition=='U':arrays=separate_chains(*arrays,chains,gap)
        audit=audit_msa(arrays[0],arrays[2],chains,gap,condition)
        assert audit['rows']<=64
        np.savez_compressed(out/f'msa_features_{calls}.npz',msa=arrays[0],deletion=arrays[1],paired=arrays[2],chain_ids=chains)
        (out/f'feature_audit_{calls}.json').write_text(json.dumps(audit,indent=2)+'\n')
        calls+=1
        return tuple(torch.as_tensor(a,dtype=t.dtype) for a,t in zip(arrays,[msa,deletion,paired]))
    features.construct_paired_msa=controlled
    native_trainer=boltz_main.Trainer
    def cpu_trainer(*args,**kwargs):
        assert kwargs['accelerator']=='cpu'
        # The native CLI selects BF16 for Boltz-2 on every accelerator. Use
        # explicit FP32 on the AVX2 CPU nodes and record this backend difference.
        kwargs['precision']='32-true'
        return native_trainer(*args,**kwargs)
    boltz_main.Trainer=cpu_trainer
    start=time.monotonic();success=False
    try:
        boltz_main.cli.main(args=sys.argv[1:],standalone_mode=False)
        assert calls>0,'No actual features were audited'
        success=True
    finally:
        metrics=dict(success=success,accelerator='cpu',precision='32-true',threads=threads,
            seconds=time.monotonic()-start,feature_audits=calls,
            peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            slurm_job_id=os.environ['SLURM_JOB_ID'],
            entry_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        (out/'execution_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')


if __name__=='__main__':main()
