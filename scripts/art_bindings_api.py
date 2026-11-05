"""Credential-safe, resumable NVIDIA transport and multi-model verification."""
from collections import Counter
import csv
import io
import json
import math
import os
from pathlib import Path
import shutil
import threading
import time
import urllib.error
import urllib.request

import numpy as np
from Bio.PDB import MMCIFParser,NeighborSearch
from Bio.SeqUtils import seq1

from art_bindings_common import ROOT,OUT,utc,sha,write,table,payload_hash,SafeRedirect

_post_lock=threading.Lock()
_last_post=0.0


def archive_gateway_failures(records,max_retries=2):
    """Explicitly retry terminal gateway responses, preserving unknown outcomes.

    Socket failures and interrupted submissions with no ID are deliberately
    excluded. A gateway retry is a new invocation; it does not recover a lost
    invocation. Only returned, sequence-verified models enter the analysis.
    The caller must hold the launcher lock.
    """
    revisions=[]
    for record in records:
        folder=OUT/'runs'/record['request_hash'];path=folder/'state.json'
        if not path.exists() or (folder/'response.json').exists():continue
        state=json.loads(path.read_text())
        if state.get('status')!='failed':continue
        errors=[json.loads(p.read_text()) for p in sorted(folder.glob('http_error_*.json'))]
        if not errors or errors[-1]['status'] not in [502,503,504]:continue
        remote_status={k.lower():v for k,v in errors[-1].get('headers',{}).items()}.get('nvcf-status','').lower()
        if state.get('request_id') and remote_status!='errored':continue
        attempts=folder/'attempts';prior=list(attempts.glob('[0-9][0-9][0-9]')) if attempts.exists() else []
        if len(prior)>=max_retries:continue
        destination=attempts/f'{len(prior)+1:03d}';destination.mkdir(parents=True)
        for source in [path,*folder.glob('http_*.json')]:shutil.move(str(source),destination/source.name)
        revision=dict(request_hash=record['request_hash'],archived_attempt=str(destination.relative_to(OUT)),
            original_http_status=errors[-1]['status'],prior_request_id=state.get('request_id'),
            prior_remote_outcome='terminal_errored_header' if remote_status=='errored' else 'unknown_no_request_ID_or_result_returned',
            action='new_invocation_of_identical_frozen_payload',updated_utc=utc(),slurm_job_id=os.environ['SLURM_JOB_ID'])
        write(destination/'retry_decision.json',revision);revisions.append(revision)
    return revisions


def normalize_response(payload,body,folder):
    """Check each model against the exact input, retaining per-model outputs."""
    count=payload['diffusion_samples'];structures=body.get('structures',[])
    if len(structures)!=count:raise ValueError(f'Expected {count} structures, received {len(structures)}')
    models=[]
    def score(name,index):
        values=body.get(name)
        if values is None:return None
        if not isinstance(values,list) or len(values)!=count:raise ValueError('Unexpected score shape: '+name)
        value=values[index]
        if value is None:return None
        if len(payload['polymers'])==1 and isinstance(value,(int,float)) and not math.isfinite(value):return None
        if not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('Nonfinite or nonnumeric score: '+name)
        return float(value)
    for index,record in enumerate(structures):
        text=record.get('structure');assert isinstance(text,str) and '_atom_site.' in text
        path=folder/f'model_{index}.cif';path.write_text(text)
        structure=MMCIFParser(QUIET=True,auth_chains=False,auth_residues=False).get_structure('prediction',io.StringIO(text))[0]
        assert len(list(structure.get_chains()))==len(payload['polymers'])
        chains=[];coordinates={};plddt={};chain_summary=[]
        for chain_index,polymer in enumerate(payload['polymers']):
            chain=structure[polymer['id']];residues=list(chain.get_residues())
            sequence=''.join(seq1(res.resname) for res in residues)
            if sequence!=polymer['sequence']:raise ValueError('Output protein sequence differs from submitted input')
            if [r.id[1] for r in residues]!=list(range(1,len(residues)+1)):
                raise ValueError('Unexpected residue numbering in returned model')
            ca=np.asarray([r['CA'].coord for r in residues]);confidence=np.asarray([r['CA'].bfactor for r in residues])
            if confidence.max(initial=0)<=1:confidence=confidence*100
            assert np.isfinite(ca).all() and np.isfinite(confidence).all()
            assert confidence.min(initial=0)>=0 and confidence.max(initial=0)<=100.001
            coordinates['chain'+str(chain_index)+'_CA']=ca
            coordinates['chain'+str(chain_index)+'_plddt']=confidence
            chain_summary.append(dict(id=polymer['id'],length=len(residues),mean_CA_plddt=float(confidence.mean())))
            chains.append(chain);plddt[polymer['id']]=confidence
        np.savez_compressed(folder/f'model_{index}.coordinates.npz',**coordinates)
        contacts={}
        if len(chains)==2:
            a,b=chains
            second_atoms=[atom for atom in b.get_atoms() if atom.element not in ['H','D']]
            search=NeighborSearch(second_atoms)
            for atom in a.get_atoms():
                if atom.element in ['H','D']:continue
                for other in search.search(atom.coord,8.0,level='A'):
                    pair=(atom.get_parent().id[1]-1,other.get_parent().id[1]-1)
                    d=float(np.linalg.norm(atom.coord-other.coord))
                    contacts[pair]=min(contacts.get(pair,float('inf')),d)
            table(folder/f'model_{index}.contacts.csv',[
                dict(chain0_residue0=i,chain1_residue0=j,min_heavy_atom_distance_A=distance)
                for (i,j),distance in sorted(contacts.items())],
                ['chain0_residue0','chain1_residue0','min_heavy_atom_distance_A'])
        scores={name:score(name,index) for name in ['protein_iptm_scores','iptm_scores','ptm_scores','confidence_scores','complex_plddt_scores','complex_iplddt_scores','complex_pde_scores','complex_ipde_scores']}
        if len(chains)==2:
            assert scores['protein_iptm_scores'] is not None and 0<=scores['protein_iptm_scores']<=1
        models.append(dict(model_index=index,structure_file=path.name,structure_sha256=sha(path),
            chain_summary=chain_summary,scores=scores,
            contacts_by_cutoff={str(cutoff):sum(d<=cutoff for d in contacts.values()) for cutoff in [4,5,8]},
            coordinates_file=f'model_{index}.coordinates.npz',contacts_file=f'model_{index}.contacts.csv' if len(chains)==2 else None))
    normalized=dict(verified=True,models=models,verified_sequences=count*len(payload['polymers']),
        server_metrics=body.get('metrics',{}),full_PAE_returned=body.get('pae') is not None,
        structures_have_unique_hashes=len({m['structure_sha256'] for m in models})==count)
    write(folder/'normalized.json',normalized)
    return normalized


def run_request(record,key,config):
    """Resume known request IDs; never blindly repeat an uncertain submission."""
    global _last_post
    folder=OUT/'runs'/record['request_hash'];folder.mkdir(parents=True,exist_ok=True)
    request_path=OUT/record['request_path'];payload=json.loads(request_path.read_text())
    assert payload_hash(payload)==record['request_hash']
    state_path=folder/'state.json';response_path=folder/'response.json'
    previous=json.loads(state_path.read_text()) if state_path.exists() else {}
    if response_path.exists():
        if (folder/'normalized.json').exists():
            normalized=json.loads((folder/'normalized.json').read_text())
            assert normalized['verified'] and len(normalized['models'])==payload['diffusion_samples']
            for model in normalized['models']:assert sha(folder/model['structure_file'])==model['structure_sha256']
        else:normalized=normalize_response(payload,json.loads(response_path.read_text()),folder)
        state=dict(previous,status='complete',completed_models=len(normalized['models']))
        write(state_path,state);return state
    if previous.get('status')=='submitting' and not previous.get('request_id'):
        raise RuntimeError('Uncertain prior submission; remote request reconciliation required')
    if previous.get('status')=='submission_uncertain' and not previous.get('request_id'):
        raise RuntimeError('Unknown prior request ID; refusing a duplicate POST')
    if previous.get('status')=='failed' and not previous.get('request_id'):
        raise RuntimeError('Previously failed request retained; inspect before explicit retry')
    request_id=previous.get('request_id');start=time.monotonic()
    state=dict(request_hash=record['request_hash'],status='running' if request_id else 'submitting',
        request_id=request_id,started_utc=previous.get('started_utc') or utc(),slurm_job_id=os.environ['SLURM_JOB_ID'],
        expected_models=payload['diffusion_samples'],stage=record['stages'],cost_precheck=False)
    write(state_path,state)
    # This hosted endpoint returns terminal 504/errored when a short queue
    # window expires. Allow the documented maximum wait before polling a
    # returned pending ID. The frozen scientific payload is unchanged.
    poll_seconds=int(config.get('poll_seconds',300))
    state['poll_seconds']=poll_seconds
    headers={'Authorization':'Bearer '+key,'Content-Type':'application/json','Accept':'application/json','NVCF-POLL-SECONDS':str(poll_seconds)}
    opener=urllib.request.build_opener(SafeRedirect())
    url='https://health.api.nvidia.com/v1/status/'+request_id if request_id else config['endpoint']
    method='GET' if request_id else 'POST';data=None if request_id else json.dumps(payload).encode()
    errors429=0;gateway_errors=0
    try:
        if method=='POST':
            with _post_lock:
                remaining=2-(time.monotonic()-_last_post)
                if remaining>0:time.sleep(remaining)
                _last_post=time.monotonic()
        for poll in range(120):
            req=urllib.request.Request(url,data=data,headers=headers,method=method)
            try:
                with opener.open(req,timeout=max(180,poll_seconds+60)) as response:
                    status=response.status
                    request_id=response.headers.get('NVCF-REQID') or request_id
                    if request_id:
                        state.update(status='running',request_id=request_id);write(state_path,state)
                    raw=response.read(100_000_001)
                    if len(raw)>100_000_000:raise ValueError('API response exceeds 100-MB limit')
                    raw=raw.replace(key.encode(),b'[REDACTED]')
                    body=json.loads(raw) if raw else {}
                    request_id=request_id or body.get('requestId') or body.get('request_id')
                    metadata={k:v for k,v in response.headers.items() if k.lower() in ['nvcf-reqid','nvcf-status','content-type','date','x-model-version']}
                    write(folder/f'http_{poll:03d}.json',dict(status=status,headers=metadata,received_utc=utc()))
            except urllib.error.HTTPError as error:
                metadata={k:v.replace(key,'[REDACTED]') for k,v in error.headers.items()
                    if k.lower() in ['nvcf-reqid','nvcf-status','content-type','date','retry-after','x-model-version']}
                request_id=error.headers.get('NVCF-REQID') or request_id
                state.update(request_id=request_id,last_http_status=error.code,
                    remote_status=error.headers.get('NVCF-STATUS'))
                write(state_path,state)
                message=error.read(16000).decode(errors='replace').replace(key,'[REDACTED]')
                write(folder/f'http_error_{poll:03d}.json',dict(status=error.code,method=method,headers=metadata,message=message,received_utc=utc()))
                if error.code==429 and errors429<8:
                    errors429+=1;time.sleep(min(60,5*2**(errors429-1)));continue
                if method=='GET' and error.code==404 and url.startswith('https://health.api.nvidia.com/'):
                    url='https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/'+request_id;continue
                if error.code in [502,503,504] and request_id and gateway_errors<5 and error.headers.get('NVCF-STATUS','').lower()!='errored':
                    # Retry only a GET for the same accepted invocation. Never
                    # issue a second POST to recover an in-flight request.
                    gateway_errors+=1;method='GET';data=None
                    url='https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/'+request_id
                    state.update(status='running');write(state_path,state)
                    time.sleep(min(30,5*gateway_errors));continue
                raise RuntimeError('HTTP '+str(error.code)+' from NVIDIA; sanitized response preserved') from None
            if status==200:
                write(response_path,body)
                normalized=normalize_response(payload,body,folder)
                state.update(status='complete',request_id=request_id,completed_models=len(normalized['models']),
                    seconds=time.monotonic()-start,completed_utc=utc(),rate_limit_responses=errors429)
                write(state_path,state)
                return state
            if status!=202 or not request_id:raise RuntimeError('Unexpected response without a resumable request identifier')
            state.update(status='running',request_id=request_id);write(state_path,state)
            method='GET';data=None;url='https://health.api.nvidia.com/v1/status/'+request_id;time.sleep(5)
        raise TimeoutError('Polling limit reached; resume the preserved request ID')
    except Exception as error:
        state.update(status='failed' if isinstance(error,(RuntimeError,ValueError,AssertionError)) else 'submission_uncertain',
            request_id=request_id,error=str(error).replace(key,'[REDACTED]'),seconds=time.monotonic()-start)
        write(state_path,state)
        raise RuntimeError(state['error']) from None
