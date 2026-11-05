"""Fetch public Boltz assets and freeze the isolated environment on a CPU node."""
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import socket
import urllib.request


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Use Slurm; no login-node computation.')
    root=Path(__file__).resolve().parents[1]
    out=root/'data/ART binding/env';out.mkdir(parents=True,exist_ok=True)
    cache=root/'.tools/boltz_cache';cache.mkdir(exist_ok=True)
    import boltz.main as boltz_main
    import torch
    assert importlib.metadata.version('boltz')=='2.2.1'
    # Make the upstream downloader's cache writes atomic so an interrupted
    # transfer cannot masquerade as a complete checkpoint on the next run.
    original=urllib.request.urlretrieve
    requests=[]
    def atomic_retrieve(url,filename=None,reporthook=None,data=None):
        destination=Path(filename);part=destination.with_suffix(destination.suffix+'.part')
        print('Fetching public model asset: '+url,flush=True)
        total=0
        with urllib.request.urlopen(url,data=data,timeout=120) as response,part.open('wb') as output:
            while True:
                chunk=response.read(8*1024*1024)
                if not chunk:break
                total+=len(chunk)
                if total>8_000_000_000:raise ValueError('Per-asset transfer limit exceeded')
                output.write(chunk)
            expected=response.headers.get('Content-Length')
            if expected:assert total==int(expected)
            requests.append(dict(url=url,response_url=response.url,bytes=total,etag=response.headers.get('ETag')))
        part.replace(destination)
        return str(destination),response.headers
    urllib.request.urlretrieve=atomic_retrieve
    try:boltz_main.download_boltz2(cache)
    finally:urllib.request.urlretrieve=original
    assets={}
    for path in sorted(cache.glob('*')):
        if path.is_file() and not path.name.endswith('.part'):
            digest=hashlib.sha256()
            with path.open('rb') as handle:
                for block in iter(lambda:handle.read(16*1024*1024),b''):digest.update(block)
            assets[path.name]=dict(bytes=path.stat().st_size,sha256=digest.hexdigest())
    for name in ['boltz2_conf.ckpt','boltz2_aff.ckpt','mols.tar']:assert name in assets
    result=dict(passed=True,slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
        boltz_version=importlib.metadata.version('boltz'),torch_version=torch.__version__,
        torch_CUDA_version=torch.version.cuda,API_key_required=False,
        optional_custom_CUDA_kernels=False,assets=assets,requests=requests,
        installed_main_source_sha256=hashlib.sha256(Path(boltz_main.__file__).read_bytes()).hexdigest())
    (out/'boltz_provenance.json').write_text(json.dumps(result,indent=2)+'\n')
    for module in ['main.py','data/parse/schema.py','data/parse/csv.py','data/feature/featurizer.py','data/feature/featurizerv2.py']:
        source=Path(boltz_main.__file__).parent/module
        if source.is_file():
            target=out/'source'/module;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
