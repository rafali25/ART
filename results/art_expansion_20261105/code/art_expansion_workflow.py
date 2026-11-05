"""Checkpointed public ART acquisition; execution requires an actual Slurm job."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import itertools
import json
import os
from pathlib import Path
import shutil
import socket
import sqlite3
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from Bio import Align, AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
import pyrodigal_gv
import zstandard

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/art_expansion_20261105'
OLD = ROOT/'results/art_sample_expansion_20261104'
CPUS = min(8, int(os.environ.get('SLURM_CPUS_PER_TASK','1')))
PROFILES = {
 'ART_RT':ROOT/'results/art_discovery/profiles/ART_RT.hmm',
 'Type_I_partner':ROOT/'results/art_discovery/profiles/Type_I_partner.hmm',
 'Type_II_partner':ROOT/'results/art_discovery/profiles/Type_II_partner.hmm',
 'TypeI_two_seed':ROOT/'results/art_type1_public/profiles/TypeI_two_seed_exploratory.hmm',
 'RVT_1':ROOT/'.tools/art-reconstruction/PF00078.hmm',
}
DATASETS = {
 'HPGC':dict(url='https://ndownloader.figshare.com/files/59221811',bytes=560950249,
    md5='03c4f8ae0ce599af7cc256c9f53a56c1',file='HPGC.7230.genomes.final.fasta.gz'),
 'JPGD':dict(url='https://zenodo.org/api/records/14974448/files/Jumbo-phage.all.fasta.gz/content',
    bytes=871869424,md5='fb500e77a3d77dd98655fdccf71e0353',file='Jumbo-phage.all.fasta.gz'),
}


def utc(): return datetime.now(timezone.utc).isoformat()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,default=str)+'\n');temp.replace(path)


def read_tsv(path):
    if not Path(path).exists():return []
    with Path(path).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def table(path,rows,fields=None):
    rows=list(rows);fields=fields or list(dict.fromkeys(k for r in rows for k in r)) or ['status']
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)


def checkpoint(n,title,summary,artifacts=()):
    folder=OUT/'checkpoints';folder.mkdir(exist_ok=True)
    if (folder/f'{n:02d}.json').exists():
        history=folder/'history';history.mkdir(exist_ok=True)
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
        for ext in ['json','md']:shutil.copy2(folder/f'{n:02d}.{ext}',history/f'{n:02d}_{stamp}.{ext}')
    obj=dict(step=n,title=title,status='complete',completed_utc=utc(),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
        summary=summary,artifacts={str(Path(p).relative_to(ROOT)):sha(p) for p in artifacts if Path(p).is_file()},
        script_sha256=sha(__file__))
    write_json(folder/f'{n:02d}.json',obj)
    lines=[f'# Checkpoint {n}: {title}','',f'Completed: {obj["completed_utc"]}',
        f'Slurm job: {obj["slurm_job_id"]}; compute host: {obj["host"]}','','```json',json.dumps(summary,indent=2),'```','']
    (folder/f'{n:02d}.md').write_text('\n'.join(lines))
    print(f'CHECKPOINT {n} COMPLETE: '+json.dumps(summary),flush=True)


def run(args,log,stdout=None):
    args=list(map(str,args));log=Path(log);log.parent.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();code=None
    try:
        with log.open('w') as err:
            if stdout:
                with Path(stdout).open('w') as dest:r=subprocess.run(args,stdout=dest,stderr=err,check=True)
            else:r=subprocess.run(args,stdout=err,stderr=subprocess.STDOUT,check=True)
        code=r.returncode
    except subprocess.CalledProcessError as e:
        code=e.returncode;raise
    finally:
        with (OUT/'commands.jsonl').open('a') as f:f.write(json.dumps(dict(args=args,log=str(log),stdout=str(stdout),returncode=code,seconds=time.monotonic()-start,job=os.environ['SLURM_JOB_ID']))+'\n')


def fetch(url,path,limit,expected_md5=None,expected_bytes=None,headers=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if expected_bytes is None or path.stat().st_size==expected_bytes:
            if expected_md5 is None or file_md5(path)==expected_md5:
                meta=path.with_suffix(path.suffix+'.request.json')
                if not meta.exists():write_json(meta,dict(url=url,bytes=path.stat().st_size,sha256=sha(path),md5=file_md5(path),
                    verified_cache=True,HTTP_headers_not_available=True,verified_utc=utc()))
                return path
        raise ValueError('Cached source validation failed: '+str(path))
    part=path.with_suffix(path.suffix+'.part')
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'ART-public-acquisition/1.0',**(headers or {})})
            total=0
            with urllib.request.urlopen(req,timeout=90) as response,part.open('wb') as out:
                info=dict(status=response.status,headers=dict(response.headers),response_url=response.url)
                while True:
                    block=response.read(min(2*1024*1024,limit-total+1))
                    if not block:break
                    total+=len(block)
                    if total>limit:raise ValueError('Byte cap exceeded')
                    out.write(block)
            if expected_bytes is not None and total!=expected_bytes:raise ValueError('Source byte count changed')
            if expected_md5 is not None and file_md5(part)!=expected_md5:raise ValueError('Published MD5 differs')
            part.replace(path)
            write_json(path.with_suffix(path.suffix+'.request.json'),dict(url=url,retrieved_utc=utc(),bytes=total,sha256=sha(path),md5=file_md5(path),**info))
            return path
        except Exception:
            if attempt==2:raise
    raise AssertionError('unreachable')


def file_md5(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def review():
    verify=json.loads((OLD/'verification.json').read_text());assert verify['passed']
    reviewed=read_tsv(OLD/'new_sample_candidates.tsv');assert len(reviewed)==4
    checked=0
    catalog=ROOT/'results/art_type1_catalog'
    old_summary=json.loads((catalog/'summary.json').read_text())
    for base,hashes in [(catalog,json.loads((catalog/'output_sha256.json').read_text())),(ROOT,old_summary['input_manifest_sha256'])]:
        for name,expected in hashes.items():assert sha(base/name.replace('\\','/'))==expected;checked+=1
    for sub in ['deferred_acquisition/processed','diverse_partner_search','jug_screen']:
        v=json.loads((OLD/sub/'verification.json').read_text())
        assert v.get('passed',v.get('source_hashes_verified',False))
    table(OUT/'reviewed_existing_additions.tsv',reviewed)
    write_json(OUT/'cluster_execution.json',dict(policy_file='ciai cluster.docx',policy_sha256=sha(ROOT/'ciai cluster.docx'),
        login_work='file preparation and Slurm administration only',execution='Slurm CPU batch',
        partition='cscc-cpu-p',qos='cscc-cpu-qos',cpus=CPUS,job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),
        no_GPU_requested=True,no_salloc=True,no_exclusive=True,no_idle_allocation=True))
    checkpoint(1,'Review existing pairs and freeze evidence tiers',dict(reviewed_pairs=4,
        passing_complete_typeI=1,passing_complete_typeII=2,RT_coverage_review_only=1,
        original_catalog_hashes_verified=checked,accepted_new_arrays=0),[OUT/'reviewed_existing_additions.tsv',OUT/'cluster_execution.json'])


def download_catalogs():
    for name,entry in DATASETS.items():
        print('Downloading '+name,flush=True)
        path=fetch(entry['url'],OUT/'downloads'/entry['file'],entry['bytes']+1,entry['md5'],entry['bytes'])
        print(name+' verified: '+str(path.stat().st_size)+' bytes',flush=True)


def screen_sources(sources,out):
    """Run the existing two-route pipeline with native executables and frozen models."""
    import importlib.util
    import art_type1_pipeline as pipeline
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    table(out/'acquired_source_records.tsv',sources)
    if not sources:
        write_json(out/'processing_summary.json',dict(source_contigs=0,primary_typeI_pairs=0));return
    driverpath=OLD/'deferred_acquisition/screen_deferred.py'
    spec=importlib.util.spec_from_file_location('art_native_recorded',driverpath)
    native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)
    cfg=json.loads((ROOT/'config/art_type1_mgnify.json').read_text());cfg.update(output=out.relative_to(ROOT).as_posix(),threads=min(4,CPUS))
    write_json(out/'run_config.json',cfg)
    pipeline.Runner=native.NativeRunner
    pipeline.profile_setup=lambda out,cfg,runner:PROFILES['TypeI_two_seed']
    original=subprocess.run
    def native_mafft(args,**kwargs):
        if args and str(args[0])=='cmd.exe':
            args=list(map(str,args));dest=Path(args[args.index('--out')+1])
            with dest.open('w') as f:r=original([str(ROOT/'.venv/bin/mafft'),'--auto','--thread','1','--inputorder',args[-1]],stdout=f,stderr=kwargs.get('stdout'),check=True)
            return r
        return original(args,**kwargs)
    subprocess.run=native_mafft
    try:pipeline.process(cfg,out,include_control=False)
    finally:subprocess.run=original
    write_json(out/'input_sha256.json',{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),driverpath,ROOT/'scripts/art_type1_pipeline.py',*PROFILES.values()]})


def expand_leads():
    from art_expansion_leads import archive_acquire,mgnify_occurrences
    archive,a=archive_acquire();mgnify,m=mgnify_occurrences()
    sources=list({r['sequence_sha256']:r for r in archive+mgnify}.values())
    screen_sources(sources,OUT/'step2/processed')
    p=OUT/'step2/processed'
    stats=json.loads((p/'processing_summary.json').read_text())
    checkpoint(2,'Expand and screen archive and MGnify occurrences',dict(archive=a,mgnify=m,screen=stats),
        [OUT/'step2/archive/frozen_leads.tsv',OUT/'step2/archive/lead_outcomes.tsv',OUT/'step2/mgnify/occurrence_manifest.tsv',p/'processing_summary.json',p/'typeI_pairs.tsv'])


def screen_catalog(name):
    from art_expansion_bulk import bulk_screen
    summary=bulk_screen(name)
    if name=='JPGD':
        from art_expansion_extra import select_uhgv
        summary=dict(summary,additional_stratum=select_uhgv())
    n=3 if name=='HPGC' else 4
    checkpoint(n,'Screen full '+name+' after exact DNA deduplication',summary,
        [OUT/name/'summary.json',OUT/name/'genome_manifest.tsv',OUT/name/'pairs.tsv',OUT/name/'profile_hits.tsv'])


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run this workload using sbatch/srun; login-node execution is prohibited.')
    OUT.mkdir(parents=True,exist_ok=True)
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','review','download','leads','HPGC','JPGD','feasibility'])
    args=parser.parse_args()
    if args.stage in ('prepare','review'):review()
    if args.stage in ('prepare','download'):download_catalogs()
    if args.stage=='leads':expand_leads()
    if args.stage in ('HPGC','JPGD'):screen_catalog(args.stage)
    if args.stage=='feasibility':
        from art_expansion_feasibility import main as assess
        assess()


if __name__=='__main__':main()
