"""Read the user-provided credential without emitting or persisting its value."""
import json
import os
from pathlib import Path
import re


def credential(root):
    text=(Path(root)/'boltz.md').read_text().strip()
    found=re.findall(r'(?:nvapi-[A-Za-z0-9_-]{20,}|sk_bc_(?:admin|ws)_(?:live|test)_[A-Za-z0-9_-]{16,})',text)
    found=list(dict.fromkeys(found))
    if len(found)!=1:raise ValueError('Expected exactly one recognized Boltz or NVIDIA API credential in boltz.md; value withheld.')
    key=found[0]
    provider='nvidia_nim' if key.startswith('nvapi-') else 'boltz_platform'
    mode='test' if provider=='boltz_platform' and '_test_' in key else 'live'
    return key,provider,mode


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm for this preparation.')
    root=Path(__file__).resolve().parents[1]
    key,provider,mode=credential(root)
    out=root/'data/ART binding/api_pilot_20261105';out.mkdir(parents=True,exist_ok=True)
    info=dict(provider=provider,credential_mode=mode,credential_value_logged=False,cost_checks_requested=False,
        slurm_job_id=os.environ['SLURM_JOB_ID'])
    (out/'provider.json').write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps(info,indent=2),flush=True)


if __name__=='__main__':main()
