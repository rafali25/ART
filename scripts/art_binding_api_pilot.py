"""Use the authorized NVIDIA Boltz-2 key for a bounded remote development run."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib
import io
import json
import os
from pathlib import Path
import socket
import time
import urllib.error
import urllib.parse
import urllib.request

from Bio.PDB import MMCIFParser,NeighborSearch

from art_binding_api_credentials import credential

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/ART binding/api_pilot_20261105'
ENDPOINT='https://health.api.nvidia.com/v1/biology/mit/boltz2/predict'


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


class SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        new=super().redirect_request(req,fp,code,msg,headers,newurl)
        if new and urllib.parse.urlsplit(req.full_url).netloc!=urllib.parse.urlsplit(newurl).netloc:
            new.remove_header('Authorization')
        return new


def payload_for(cell):
    path=ROOT/cell['input_path']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==cell['input_sha256']
    original=json.loads(path.read_text());assert cell['condition']=='S'
    polymers=[]
    for entity in original['sequences']:
        p=entity['protein'];assert p['msa']=='empty'
        sequence=p['sequence']
        polymers.append(dict(id=p['id'],molecule_type='protein',sequence=sequence,
            msa={'query_only':{'a3m':dict(format='a3m',alignment='>query\n'+sequence+'\n',rank=0)}}))
    return dict(polymers=polymers,ligands=[],constraints=[],recycling_steps=3,sampling_steps=200,
        diffusion_samples=1,step_scale=1.5,without_potentials=True,output_format='mmcif',concatenate_msas=False)


def collect_result(cell,body,folder):
    if not isinstance(body,dict) or not body.get('structures'):
        raise ValueError('API response contains no predicted structure; inspect the saved response.')
    assert len(body['structures'])==1,'One diffusion sample was requested'
    text=body['structures'][0]['structure'];assert isinstance(text,str) and '_atom_site.' in text
    path=folder/'model_0.cif';path.write_text(text)
    structure=MMCIFParser(QUIET=True,auth_chains=False,auth_residues=False).get_structure('prediction',io.StringIO(text))[0]
    chains=list(structure.get_chains());assert len(chains)==2
    lengths=[len(list(c.get_residues())) for c in chains]
    assert lengths==[cell['RT_length_aa'],cell['partner_length_aa']],str(lengths)
    atoms=[a for a in structure.get_atoms() if a.element not in ['H','D']]
    contacts=set()
    for a,b in NeighborSearch(atoms).search_all(5.0,level='A'):
        ra,rb=a.get_parent(),b.get_parent()
        if ra.get_parent().id!=rb.get_parent().id:
            if ra.get_parent().id!=chains[0].id:ra,rb=rb,ra
            contacts.add((ra.id[1],rb.id[1]))
    score=body.get('protein_iptm_scores',body.get('iptm_scores'))
    assert score and len(score)==1 and 0<=score[0]<=1
    result=dict(cell_id=cell['cell_id'],RT_sample_id=cell['RT_sample_id'],partner_sample_id=cell['partner_sample_id'],
        native_pair=cell['native_pair'],provider='NVIDIA_hosted_Boltz2',condition='S_explicit_query_only_MSA',
        protein_iptm=score[0],iptm=body.get('iptm_scores'),ptm=body.get('ptm_scores'),
        confidence=body.get('confidence_scores'),contacts_5A_residue_pairs=len(contacts),
        structure_file=str(path.relative_to(ROOT)),structure_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        chain_lengths=lengths,server_runtime_metrics=body.get('metrics',{}),
        exact_seed_control=False,processed_features_inspectable=False,full_PAE_available=body.get('pae') is not None,
        confirmatory=False,experimental_binding_validation=False)
    write(folder/'result.json',result)
    return result


def run_cell(cell,key):
    folder=OUT/cell['cell_id'];folder.mkdir(exist_ok=True)
    state_path=folder/'state.json';response_path=folder/'response.json'
    if response_path.exists():
        result=collect_result(cell,json.loads(response_path.read_text()),folder)
        return result
    previous=json.loads(state_path.read_text()) if state_path.exists() else {}
    if previous.get('status') in ['submitting','submission_uncertain'] and not previous.get('request_id'):
        raise RuntimeError('Uncertain prior submission; refusing an automatic duplicate for '+cell['cell_id'])
    payload=payload_for(cell);write(folder/'request.json',payload)
    state=dict(cell_id=cell['cell_id'],status='submitting',endpoint=ENDPOINT,
        started_utc=datetime.now(timezone.utc).isoformat(),slurm_job_id=os.environ['SLURM_JOB_ID'],
        request_sha256=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest(),cost_precheck=False,
        seed_parameter_available=False,submitted_seed=None,registered_local_seed=cell['seed'])
    request_id=previous.get('request_id');start=time.monotonic();write(state_path,state)
    headers={'Authorization':'Bearer '+key,'Content-Type':'application/json','Accept':'application/json','NVCF-POLL-SECONDS':'60'}
    opener=urllib.request.build_opener(SafeRedirect())
    url='https://health.api.nvidia.com/v1/status/'+request_id if request_id else ENDPOINT
    method='GET' if request_id else 'POST'
    data=None if request_id else json.dumps(payload).encode()
    print('Submitting '+cell['cell_id']+' to NVIDIA Boltz-2',flush=True)
    try:
        for poll in range(90):
            request=urllib.request.Request(url,data=data,headers=headers,method=method)
            try:
                with opener.open(request,timeout=180) as response:
                    status=response.status;raw=response.read(50_000_001)
                    if len(raw)>50_000_000:raise ValueError('Unexpectedly large API response')
                    metadata={k:v for k,v in response.headers.items() if k.lower() in ['nvcf-reqid','content-type','date','nvcf-status','x-model-version']}
                    body=json.loads(raw) if raw else {}
                    write(folder/f'http_{poll:03d}.json',dict(status=status,headers=metadata,
                        received_utc=datetime.now(timezone.utc).isoformat(),body=body if status!=200 else {'saved_in':'response.json'}))
                    request_id=response.headers.get('NVCF-REQID') or request_id or body.get('requestId') or body.get('request_id')
            except urllib.error.HTTPError as error:
                text=error.read(16000).decode(errors='replace').replace(key,'[REDACTED]')
                write(folder/f'http_error_{poll:03d}.json',dict(status=error.code,message=text))
                if error.code==429 and poll<3:
                    time.sleep(10);continue
                if method=='GET' and error.code==404 and url.startswith('https://health.api.nvidia.com/'):
                    url='https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/'+request_id;continue
                raise RuntimeError('NVIDIA API returned HTTP '+str(error.code)+'; response saved without credentials') from None
            if status==200:
                write(response_path,body)
                result=collect_result(cell,body,folder)
                state.update(status='complete',request_id=request_id,seconds=time.monotonic()-start,
                    ended_utc=datetime.now(timezone.utc).isoformat())
                write(state_path,state);result['client_elapsed_seconds']=state['seconds'];write(folder/'result.json',result)
                print('Completed '+cell['cell_id']+' in '+str(round(state['seconds'],1))+' seconds',flush=True)
                return result
            if status!=202 or not request_id:raise RuntimeError('Unexpected asynchronous response without a request identifier')
            state.update(status='running',request_id=request_id);write(state_path,state)
            url='https://health.api.nvidia.com/v1/status/'+request_id;data=None;method='GET'
            time.sleep(5)
        raise TimeoutError('API polling limit reached; request identifier preserved for resume')
    except Exception as error:
        state.update(status='failed' if isinstance(error,RuntimeError) else 'submission_uncertain',request_id=request_id,
            seconds=time.monotonic()-start,error=str(error).replace(key,'[REDACTED]'))
        write(state_path,state)
        raise RuntimeError(state['error']) from None


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm for API requests and analysis.')
    key,provider,mode=credential(ROOT)
    assert provider=='nvidia_nim' and mode=='live','This runner requires the recognized NVIDIA hosted-inference key'
    OUT.mkdir(exist_ok=True)
    manifest=json.loads((ROOT/'data/ART binding/pilot_20261105/manifest.json').read_text())
    cells=[c for c in manifest['cells'] if c['condition']=='S'];assert len(cells)==4
    protocol=dict(provider='NVIDIA_hosted_Boltz2',endpoint=ENDPOINT,host=socket.gethostname(),
        slurm_job_id=os.environ['SLURM_JOB_ID'],cells=cells,condition='explicit_query_only_MSA',
        purpose='remote_development_pilot',cost_checks=False,credential_logged=False,
        sampling_steps=200,recycling_steps=3,diffusion_samples=1,step_scale=1.5,parallel_requests=2,
        no_templates=True,no_restraints=True,no_affinity_module=True,
        limitations=['NVIDIA endpoint does not expose a seed parameter: samples are not registered fixed-seed replicates.',
            'Remote processed features and exact checkpoint hash are not inspectable; outputs remain exploratory.',
            'Strict block-diagonal unpaired U adapter cannot be verified remotely; U is not submitted as equivalent.',
            'Keep hosted outputs separate from local GPU/CPU outputs and calibrate the chosen backend before inference.'])
    write(OUT/'protocol.json',protocol)
    results=[];failures=[]
    with ThreadPoolExecutor(max_workers=2) as executor:
        pending={executor.submit(run_cell,cell,key):cell for cell in cells}
        for future in as_completed(pending):
            try:results.append(future.result())
            except Exception as error:
                failures.append(dict(cell_id=pending[future]['cell_id'],error=str(error).replace(key,'[REDACTED]')))
            write(OUT/'summary.json',dict(status='running',completed=len(results),planned=4,failures=failures,results=results,confirmatory=False))
    summary=dict(status='complete' if len(results)==4 else 'incomplete',completed=len(results),planned=4,
        failures=failures,results=results,provider='NVIDIA_hosted_Boltz2',cost_checks=False,credential_logged=False,
        confirmatory=False,slurm_job_id=os.environ['SLURM_JOB_ID'])
    write(OUT/'summary.json',summary)
    print(json.dumps(dict(status=summary['status'],completed=len(results),planned=4,failures=failures),indent=2),flush=True)


if __name__=='__main__':main()
