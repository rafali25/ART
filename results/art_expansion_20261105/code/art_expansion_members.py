"""Bounded real MGnify Parquet member/occurrence query through a metered proxy."""
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
import urllib.request
import urllib.parse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from art_expansion_workflow import ROOT,OUT,OLD,fetch,sha,table,read_tsv,write_json,run,checkpoint,screen_sources
from art_type1_public_search import locate_peptide


def main():
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm.')
    folder=OUT/'step2/parquet';folder.mkdir(parents=True,exist_ok=True)
    # Source BioSample/study metadata helps avoid treating repeated sequencing
    # runs as independent biological samples; sample taxon is not a phage host.
    archived=read_tsv(OUT/'step2/archive/frozen_leads.tsv')
    runs=sorted({r['run'] for r in archived}|{'SRR8925777','SRR26805139','ERR8055547','ERR10896470','SRR7532897'});contexts=[];context_errors=[]
    for start in range(0,len(runs),15):
        batch=runs[start:start+15]
        query=' OR '.join('run_accession="'+r+'"' for r in batch)
        url='https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(dict(result='read_run',query=query,
            fields='run_accession,study_accession,sample_accession,scientific_name,experiment_title',format='tsv',limit=1000))
        try:
            path=fetch(url,folder/f'ENA_runs_{start:03d}.tsv',1000000)
            contexts.extend(read_tsv(path))
        except Exception as e:context_errors.append(dict(runs_json=json.dumps(batch),reason=str(e)))
    table(folder/'source_run_context.tsv',contexts);table(folder/'source_context_failures.tsv',context_errors)
    hits=read_tsv(ROOT/'results/art_type1_public/remote_hits.tsv')
    reps=sorted({r['protein_accession'] for r in hits if r['family']=='Type_I_partner' and r['protein_accession'].startswith('MGYP')})
    # Four explicit cluster representatives keep the remote benchmark narrow.
    chosen=reps[:4];numbers=[int(re.sub(r'^MGYP','',r)) for r in chosen]
    state=dict(bytes=0,limit=100_000_000,requests=[],lock=threading.Lock())
    allowed={'mgy_cluster_seqs.parquet','mgy_proteins_metadata.parquet','mgy_protein_sequences.parquet'}
    base='https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/2026_07/'
    class Proxy(BaseHTTPRequestHandler):
        def log_message(self,*args):pass
        def do_HEAD(self):self.forward(True)
        def do_GET(self):self.forward(False)
        def forward(self,head):
            name=self.path.lstrip('/').split('?')[0]
            if name not in allowed:self.send_error(404);return
            headers={}
            if self.headers.get('Range'):headers['Range']=self.headers['Range']
            try:
                request=urllib.request.Request(base+name,headers=headers,method='HEAD' if head else 'GET')
                with urllib.request.urlopen(request,timeout=30) as response:
                    length=int(response.headers.get('Content-Length','0'))
                    with state['lock']:
                        if not head and (length<=0 or state['bytes']+length>state['limit']):raise ValueError('Parquet remote-read byte budget exceeded')
                        if not head:state['bytes']+=length
                    body=b'' if head else response.read(length+1)
                    if not head and len(body)>length:raise ValueError('Response body exceeded declared length')
                    self.send_response(response.status)
                    for key in ['Content-Length','Content-Type','Content-Range','Accept-Ranges','ETag','Last-Modified']:
                        if response.headers.get(key):self.send_header(key,response.headers[key])
                    self.end_headers()
                    if not head:self.wfile.write(body)
                    with state['lock']:state['requests'].append(dict(file=name,range=headers.get('Range'),method='HEAD' if head else 'GET',bytes=len(body),status=response.status))
            except Exception as e:
                with state['lock']:state['requests'].append(dict(file=name,range=headers.get('Range'),error=str(e)))
                try:self.send_error(503,str(e))
                except Exception:pass
    server=None;connection=None;summary=dict(chosen_representatives=chosen,byte_cap=state['limit'],status='not_started')
    sources=[];mappings=[]
    try:
        deps=OUT/'dependencies';deps.mkdir(exist_ok=True)
        sys.path.insert(0,str(deps))
        try:import duckdb
        except ImportError:
            run([sys.executable,'-m','pip','install','--disable-pip-version-check','--no-cache-dir','--target',deps,'duckdb==1.4.3'],folder/'duckdb-install.log')
            import duckdb
        summary['duckdb_version']=duckdb.__version__
        connection=duckdb.connect();extensions=deps/'duckdb_extensions';extensions.mkdir(exist_ok=True)
        connection.execute("SET extension_directory = '"+str(extensions)+"'")
        connection.execute("INSTALL httpfs FROM 'https://extensions.duckdb.org'; LOAD httpfs;")
        connection.execute('SET threads=2; SET http_retries=1; SET http_timeout=25;')
        server=ThreadingHTTPServer(('127.0.0.1',0),Proxy)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        local=f'http://127.0.0.1:{server.server_port}/'
        def query(sql,name):
            timer=threading.Timer(150,connection.interrupt);timer.start()
            try:
                result=connection.execute(sql);columns=[c[0] for c in result.description]
                rows=[dict(zip(columns,r)) for r in result.fetchall()]
                table(folder/name,rows,columns);return rows
            finally:timer.cancel()
        members=query("SELECT cluster_rep,cluster_member FROM read_parquet('"+local+"mgy_cluster_seqs.parquet') WHERE cluster_rep IN ("+','.join(map(str,numbers))+") LIMIT 500",'members.tsv')
        member_ids=sorted({int(r['cluster_member']) for r in members})
        summary.update(cluster_members_returned=len(member_ids),member_results_capped=len(members)==500)
        if member_ids:
            occurrences=query("SELECT protein_id,contig_id,study_accession,assembly_accession,contig_name,contig_length,start_position,end_position,strand,truncation,biome FROM read_parquet('"+local+"mgy_proteins_metadata.parquet') WHERE protein_id IN ("+','.join(map(str,member_ids))+") LIMIT 500",'occurrences.tsv')
            summary['occurrences_returned']=len(occurrences)
            selected=sorted(occurrences,key=lambda r:(str(r['truncation']) not in ['0','00'],int(r['contig_length'])<15000,r['study_accession'],r['contig_id']))[:30]
            ids=sorted({int(r['protein_id']) for r in selected})
            peptides=query("SELECT protein_id,full_length,sequence FROM read_parquet('"+local+"mgy_protein_sequences.parquet') WHERE protein_id IN ("+','.join(map(str,ids))+")",'member_sequences.tsv') if ids else []
            peptide_map={int(r['protein_id']):r['sequence'] for r in peptides}
            known={r['accession'] for c in ['art_type1_public','art_type1_mgnify'] for r in read_tsv(ROOT/'results'/c/'source_records.tsv')}
            known.update(r['accession'] for r in read_tsv(OUT/'step2/archive/acquired_source_records.tsv'))
            seen=set();body_bytes=0;failures=[]
            for loc in selected:
                cid=f'MGYC{int(loc["contig_id"]):012d}'
                if cid in known or cid in seen or int(loc['protein_id']) not in peptide_map:continue
                if len(seen)>=10 or body_bytes>=20_000_000:break
                seen.add(cid)
                try:
                    assembly=loc['assembly_accession']
                    ap=fetch('https://www.ebi.ac.uk/metagenomics/api/v2/assemblies/'+assembly+'/analyses?page_size=100',folder/(assembly+'.json'),1000000)
                    data=json.loads(ap.read_text());analyses=data.get('items',data.get('results',[]));record=None
                    for a in analyses[:3]:
                        aid=a.get('accession') or a.get('mgya_accession')
                        url='https://www.ebi.ac.uk/metagenomics/api/v1/analyses/'+aid+'/contigs/'+urllib.parse.quote(loc['contig_name'],safe='')
                        try:
                            path=fetch(url,folder/(cid+'.fna'),min(5_000_000,20_000_000-body_bytes),headers={'Accept':'application/json'})
                            body_bytes+=path.stat().st_size;record=SeqIO.read(path,'fasta');break
                        except Exception as e:last_error=str(e)
                    if record is None:raise ValueError(last_error if analyses else 'No analysis mapping')
                    if record.id!=loc['contig_name']:raise ValueError('Contig metadata identity differs')
                    record.id=cid;record.annotations.update(molecule_type='DNA',topology='linear')
                    positions=locate_peptide(record,peptide_map[int(loc['protein_id'])],int(loc['strand']))
                    mappings.append(dict(member_protein=f'MGYP{int(loc["protein_id"]):012d}',contig=cid,maps_json=json.dumps(positions),unique=len(positions)==1))
                    if len(positions)!=1:continue
                    dest=folder/(cid+'.gb');SeqIO.write(record,dest,'genbank')
                    sources.append(dict(accession=cid,path=dest.relative_to(ROOT).as_posix(),format='genbank',file_sha256=sha(dest),
                        sequence_sha256=hashlib.sha256(str(record.seq).upper().encode()).hexdigest(),source_url=url,
                        source_version_kind='MGnify_2026_07_Parquet_member_snapshot',database_release='2026_07',role='public_cluster_member_occurrence',bytes=path.stat().st_size))
                except Exception as e:failures.append(dict(contig=cid,reason=str(e)))
            table(folder/'source_failures.tsv',failures);summary.update(new_source_contigs=len(sources),new_source_bytes=body_bytes)
        summary['status']='bounded_query_complete'
    except Exception as e:
        summary.update(status='bounded_query_unavailable_or_budget_limited',reason=type(e).__name__+': '+str(e),
            limitation='A failed or capped remote column scan is not absence of members or ART loci')
    finally:
        if connection:connection.close()
        if server:server.shutdown();server.server_close()
    table(folder/'acquired_source_records.tsv',sources);table(folder/'peptide_mappings.tsv',mappings)
    write_json(folder/'remote_read_audit.json',dict(bytes=state['bytes'],limit=state['limit'],requests=state['requests']))
    summary['remote_Parquet_bytes']=state['bytes'];write_json(folder/'summary.json',summary)
    if sources:screen_sources(sources,folder/'processed')
    prior=json.loads((OUT/'checkpoints/02.json').read_text())['summary']
    checkpoint(2,'Archive and MGnify occurrence/member expansion',dict(prior,Parquet_member_expansion=summary),
        [folder/'summary.json',folder/'remote_read_audit.json',folder/'members.tsv',folder/'occurrences.tsv',folder/'peptide_mappings.tsv'])


if __name__=='__main__':main()
