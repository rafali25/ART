"""Frozen inputs and resumable transport for exploratory ART interface studies."""
from datetime import datetime,timezone
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request
import urllib.parse

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/results_bindings'


def utc():return datetime.now(timezone.utc).isoformat()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def table(path,rows,fields=None):
    rows=list(rows);fields=fields or (list(rows[0]) if rows else ['status'])
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)


def read_table(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))


def cropped_a3m(query_aligned,homolog_aligned,start0,end0):
    """Project an aligned homolog onto a query interval, preserving insertions."""
    if len(query_aligned)!=len(homolog_aligned):raise ValueError('Alignment lengths differ')
    query_position=0;result=[]
    for q,h in zip(query_aligned,homolog_aligned):
        if q!='-':
            if start0<=query_position<end0:result.append(h)
            query_position+=1
        elif start0<query_position<end0 and h!='-':result.append(h.lower())
    if not 0<=start0<end0<=query_position:raise ValueError('Invalid query interval')
    return ''.join(result)


def msa_csv(sequences):
    text=io.StringIO();writer=csv.writer(text,lineterminator='\n');writer.writerow(['sequence','key'])
    for sequence in sequences:writer.writerow([sequence,''])
    return text.getvalue()


def polymer(chain_id,sequence,condition,msa_path=None):
    if condition=='S_query_only':
        alignment='>query\n'+sequence+'\n';format='a3m'
    elif condition=='U_requested_independent_MSA_order':
        if not msa_path:raise ValueError('U-requested requires a frozen per-chain MSA')
        alignment=Path(msa_path).read_text();format='csv'
        rows=list(csv.DictReader(io.StringIO(alignment)))
        if not rows or rows[0]['sequence']!=sequence:raise ValueError('MSA query differs from input')
        if any(r['key'] for r in rows):raise ValueError('Pairing identifiers must be blank')
        # The hosted NVIDIA schema is stricter than native Boltz's CSV reader:
        # it requires this exact column order, although the values are unchanged.
        stream=io.StringIO();writer=csv.writer(stream,lineterminator='\n');writer.writerow(['key','sequence'])
        for row in rows:writer.writerow(['',row['sequence']])
        alignment=stream.getvalue()
    else:raise ValueError('Unknown input condition')
    return dict(id=chain_id,molecule_type='protein',sequence=sequence,
        msa={'controlled_input':{format:dict(format=format,alignment=alignment,rank=0)}})


def payload_hash(payload):
    return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()


class SafeRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self,req,fp,code,msg,headers,newurl):
        result=super().redirect_request(req,fp,code,msg,headers,newurl)
        if result and urllib.parse.urlsplit(req.full_url).netloc!=urllib.parse.urlsplit(newurl).netloc:
            result.remove_header('Authorization')
        return result
