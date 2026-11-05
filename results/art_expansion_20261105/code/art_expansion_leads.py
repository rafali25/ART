"""Bounded current-Logan peptide mapping and MGnify occurrence expansion."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.request
import urllib.parse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import zstandard

from art_expansion_workflow import ROOT, OUT, OLD, sha, table, read_tsv, write_json, fetch, utc
from art_type1_public_search import locate_peptide


def records_from_prefix(raw,complete_object):
    reader=zstandard.ZstdDecompressor().stream_reader(io.BytesIO(raw),read_across_frames=True)
    chunks=[];size=0
    while True:
        chunk=reader.read(1024*1024)
        if not chunk:break
        size+=len(chunk)
        if size>300_000_000:raise ValueError('Decoded prefix exceeds 300-MB safety cap')
        chunks.append(chunk)
    text=b''.join(chunks).decode()
    if not complete_object:
        pos=text.rfind('\n>');text=text[:pos+1] if pos>=0 else ''
    return SeqIO.parse(io.StringIO(text),'fasta')


def match_records(raw,complete_object,leads,peptides):
    matches=defaultdict(list);count=0;last=''
    # Reverse-translate a short internal peptide into an exact degenerate DNA
    # regular expression. The regex is only a fast prefilter; acceptance still
    # requires locate_peptide's full exact genomic translation and uniqueness.
    codons=defaultdict(set)
    for code in [11,15]:
        for codon,amino in CodonTable.unambiguous_dna_by_id[code].forward_table.items():codons[amino].add(codon)
    probes={lead:re.compile(''.join('(?:'+'|'.join(sorted(codons[amino]))+')' for amino in str(peptides[lead].seq)[10:30])) for lead in leads}
    for rec in records_from_prefix(raw,complete_object):
        count+=1;last=rec.id;possible=set()
        for dna in [str(rec.seq).upper(),str(rec.seq.reverse_complement()).upper()]:
            possible.update(lead for lead,probe in probes.items() if probe.search(dna))
        for lead in possible:
            positions=locate_peptide(rec,str(peptides[lead].seq))
            if positions:matches[lead].append((rec,positions))
    return matches,count,last


def archive_acquire():
    folder=OUT/'step2/archive';folder.mkdir(parents=True,exist_ok=True)
    leads=[r for r in read_tsv(OLD/'deferred_leads.tsv') if r['source']=='archived_Logan_partner_sequence' and not r['retrieved_current_contig']]
    assert len(leads)==45
    peptides={r.id:r for r in SeqIO.parse(ROOT/'results/art_type1_public/archive_leads/leads_and_GNAT_controls.faa','fasta')}
    screen=read_tsv(ROOT/'results/art_type1_public/archive_leads/profile_cross_screen.tsv')
    bylead=defaultdict(list)
    for h in screen:
        if float(h['sequence_evalue'])<=1e-5 and float(h['domain_evalue'])<=1e-5 and float(h['profile_coverage'])>=.6:
            bylead[h['target']].append(h['family'])
    for r in leads:
        r['run']=r['lead_id'].split('_')[0]
        r['historical_contig_number']=int(r['lead_id'].split('_')[1])
        r['passing_frozen_profiles']=';'.join(sorted(set(bylead[r['lead_id']])))
    leads.sort(key=lambda r:(not('Type_I_partner' in r['passing_frozen_profiles']),r['historical_contig_number']>10000,
        float(r['rough_identity'])>.9,r['historical_contig_number'],float(r['rough_identity'])))
    for i,r in enumerate(leads,1):r['execution_priority']=i
    table(folder/'frozen_leads.tsv',leads)
    groups=defaultdict(list)
    for r in leads:groups[r['run']].append(r['lead_id'])
    requests=json.loads((folder/'requests.json').read_text()) if (folder/'requests.json').exists() else []
    sources={r['accession']:r for r in read_tsv(folder/'acquired_source_records.tsv')}
    results={r['lead_id']:r for r in read_tsv(folder/'lead_outcomes.tsv')}
    total=max(sum(int(r.get('body_bytes',0)) for r in requests),sum(p.stat().st_size for p in folder.glob('*.prefix.zst')))
    def prefetch(item):
        run,wanted=item;path=folder/(run+'.prefix.zst');mp=folder/(run+'.prefix.json')
        old=path.read_bytes() if path.exists() else b''
        meta=json.loads(mp.read_text()) if mp.exists() else {}
        required=wanted-len(old);url=f'https://s3.amazonaws.com/logan-pub/c/{run}/{run}.contigs.fa.zst'
        headers={'Range':f'bytes={len(old)}-{wanted-1}','User-Agent':'ART-mapped-locus-study/1.0'}
        if old:headers['If-Match']=meta['etag']
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=45) as response:
                status=response.status;rh=dict(response.headers);cr=response.headers.get('Content-Range','')
                if old and (status!=206 or not cr.startswith(f'bytes {len(old)}-')):raise ValueError('Extension range not honored')
                if status==206 and not cr.startswith(f'bytes {len(old)}-'):raise ValueError('Incorrect offset')
                extra=response.read(required+1)
                if len(extra)>required:raise ValueError('Byte cap exceeded')
                etag=response.headers.get('ETag','')
            raw=old+extra;size=int(cr.rsplit('/',1)[1]) if '/' in cr and cr.rsplit('/',1)[1].isdigit() else len(raw) if status==200 else None
            path.write_bytes(raw)
            write_json(mp,dict(run=run,url=url,etag=etag,bytes=len(raw),complete_object=size==len(raw),object_bytes=size,retrieved_utc=utc(),sha256=hashlib.sha256(raw).hexdigest()))
            print(f'Fetched {run}: {len(raw)} prefix bytes',flush=True)
            return dict(run=run,range_requested=headers['Range'],status=status,response_headers=rh,body_bytes=len(extra),time=utc())
        except Exception as e:
            print(f'Range attempt failed {run}: {e}',flush=True)
            return dict(run=run,range_target=wanted,status='failed',error=str(e),time=utc())
    # All runs get an initial attempt; later ranges improve context without
    # interpreting an incomplete stream as absence. Shared runs are downloaded once.
    for target_bytes in [1_000_000,3_000_000,6_000_000,10_000_000]:
        tasks=[];reserved=total
        for run,ids in groups.items():
            if all(results.get(lead,{}).get('status')=='exact_unique_peptide_mapping' for lead in ids):continue
            path=folder/(run+'.prefix.zst');mp=folder/(run+'.prefix.json')
            meta=json.loads(mp.read_text()) if mp.exists() else {}
            existing=path.stat().st_size if path.exists() else 0
            required=target_bytes-existing
            if meta.get('complete_object') or required<=0 or reserved+required>200_000_000:continue
            reserved+=required;tasks.append((run,target_bytes))
        print(f'Archive transfer tier {target_bytes}: {len(tasks)} tasks, max three simultaneous requests',flush=True)
        with ThreadPoolExecutor(max_workers=3) as executor:
            for request in executor.map(prefetch,tasks):
                requests.append(request);total+=int(request.get('body_bytes',0))
                write_json(folder/'requests.json',requests)
        for run,ids in groups.items():
            if all(results.get(lead,{}).get('status')=='exact_unique_peptide_mapping' for lead in ids):continue
            path=folder/(run+'.prefix.zst');meta_path=folder/(run+'.prefix.json')
            existing=path.read_bytes() if path.exists() else b''
            meta=json.loads(meta_path.read_text()) if meta_path.exists() else {}
            already_scanned=all(lead in results and int(results[lead].get('prefix_bytes') or 0)>=len(existing) for lead in ids)
            if meta.get('complete_object') and already_scanned:continue
            required=target_bytes-len(existing)
            if required<=0 and already_scanned:continue
            required=max(0,required)
            if total+required>200_000_000:continue
            url=f'https://s3.amazonaws.com/logan-pub/c/{run}/{run}.contigs.fa.zst'
            try:
                if not meta or len(existing)<target_bytes:
                    headers={'Range':f'bytes={len(existing)}-{target_bytes-1}','User-Agent':'ART-mapped-locus-study/1.0'}
                    if existing:headers['If-Match']=meta['etag']
                    req=urllib.request.Request(url,headers=headers)
                    with urllib.request.urlopen(req,timeout=45) as response:
                        code=response.status;rh=dict(response.headers);cr=response.headers.get('Content-Range','')
                        if existing and (code!=206 or not cr.startswith(f'bytes {len(existing)}-')):raise ValueError('Extension range not honored')
                        if not existing and code not in [200,206]:raise ValueError('Invalid HTTP status')
                        if code==206 and not cr.startswith(f'bytes {len(existing)}-'):raise ValueError('Wrong range offset')
                        extra=response.read(required+1)
                        if len(extra)>required:raise ValueError('Server exceeded requested range')
                        etag=response.headers.get('ETag','')
                    total+=len(extra);raw=existing+extra
                    object_size=int(cr.rsplit('/',1)[1]) if '/' in cr and cr.rsplit('/',1)[1].isdigit() else len(raw) if code==200 else None
                    meta=dict(run=run,url=url,etag=etag,bytes=len(raw),complete_object=object_size==len(raw),object_bytes=object_size,
                        retrieved_utc=utc(),sha256=hashlib.sha256(raw).hexdigest())
                    path.write_bytes(raw);write_json(meta_path,meta)
                    requests.append(dict(run=run,range_requested=headers['Range'],status=code,response_headers=rh,body_bytes=len(extra),time=utc()))
                else:raw=existing
                matches,n,last=match_records(raw,meta['complete_object'],ids,peptides)
                for lead in ids:
                    found=matches.get(lead,[])
                    result=dict(lead_id=lead,run=run,prefix_bytes=len(raw),complete_records_searched=n,last_complete_record=last,
                        complete_object=meta['complete_object'],matches=len(found))
                    if len(found)==1 and len(found[0][1])==1:
                        rec,maps=found[0];rec.annotations.update(molecule_type='DNA',topology='linear')
                        rec.description+='; current Logan object exact archived peptide mapping'
                        dest=folder/(rec.id+'.gb');SeqIO.write(rec,dest,'genbank')
                        sources[rec.id]=dict(accession=rec.id,path=dest.relative_to(ROOT).as_posix(),format='genbank',
                            file_sha256=sha(dest),sequence_sha256=hashlib.sha256(str(rec.seq).upper().encode()).hexdigest(),
                            source_url=url,source_version_kind='Logan_object_ETag_and_content_SHA256',
                            database_release='current_object;'+meta['etag'],role='step2_archive_public_lead',length_nt=len(rec),bytes=len(raw))
                        result.update(status='exact_unique_peptide_mapping',source_contig=rec.id,length_nt=len(rec),coordinates_json=json.dumps(maps))
                    elif found:result.update(status='multiple_peptide_mappings_requires_review',coordinates_json=json.dumps([dict(contig=r.id,maps=m) for r,m in found]))
                    else:result['status']='peptide_absent_from_complete_current_object' if meta['complete_object'] else 'unresolved_bounded_prefix'
                    results[lead]=result
                print(f'Archive {run}: prefix={len(raw)}, records={n}, mapped={sum(len(matches.get(i,[]))==1 for i in ids)}, total_bytes={total}',flush=True)
            except Exception as e:
                requests.append(dict(run=run,range_target=target_bytes,status='failed',error=type(e).__name__+': '+str(e),time=utc()))
                for lead in ids:results.setdefault(lead,dict(lead_id=lead,run=run,status='retrieval_or_scan_failed',reason=str(e)))
            table(folder/'lead_outcomes.tsv',list(results.values()))
            table(folder/'acquired_source_records.tsv',list(sources.values()))
            write_json(folder/'requests.json',requests)
    summary=dict(lead_rows=len(leads),distinct_source_runs=len(groups),source_contigs=len(sources),
        exact_mapped_leads=sum(r['status']=='exact_unique_peptide_mapping' for r in results.values()),
        statuses=dict(Counter(r['status'] for r in results.values())),successful_new_body_bytes=total,byte_cap=200000000,
        every_run_attempted=all(any(q['run']==run for q in requests) for run in groups),prefix_nonmatches_are_not_absence=True)
    write_json(folder/'summary.json',summary)
    return list(sources.values()),summary


def mgnify_occurrences():
    folder=OUT/'step2/mgnify';folder.mkdir(parents=True,exist_ok=True)
    hits=read_tsv(ROOT/'results/art_type1_public/remote_hits.tsv')+read_tsv(ROOT/'results/art_type1_mgnify/remote_hits.tsv')
    reps=sorted({h['protein_accession'] for h in hits if h['family']=='Type_I_partner' and h['protein_accession'].startswith('MGYP')})
    locations=[];failures=[]
    known={r['accession'] for cohort in ['art_type1_public','art_type1_mgnify'] for r in read_tsv(ROOT/'results'/cohort/'source_records.tsv')}
    known.update(r['accession'] for r in read_tsv(OLD/'diverse_partner_search/acquired_source_records.tsv'))
    known.update(r['accession'] for r in read_tsv(OLD/'deferred_acquisition/processed/acquired_source_records.tsv'))
    known.update(r['source_genome'] for r in read_tsv(OLD/'deferred_acquisition/processed/primary_pair_summary.tsv'))
    for rep in reps:
        try:
            path=fetch('https://www.ebi.ac.uk/metagenomics/proteins/api/v1/protein/'+rep,folder/(rep+'.json'),1000000)
            data=json.loads(path.read_text())
            for loc in data.get('study_assembly_contigs',[]):locations.append(dict(protein=rep,peptide=data['sequence'],full_length=data.get('full_length'),cluster_size=data.get('cluster_size'),**loc))
        except Exception as e:failures.append(dict(protein=rep,status='representative_metadata_failed',reason=str(e)))
    locations.sort(key=lambda r:(r['contig_accession'] in known,not bool(r['full_length']),min(int(r['start']),int(r['end']))<6000,r['study_accession'],r['protein']))
    table(folder/'occurrence_manifest.tsv',[{k:v for k,v in r.items() if k!='peptide'} for r in locations])
    sources={};maps=[];used=0;seen=set()
    for loc in locations:
        contig=loc['contig_accession']
        if contig in known or contig in seen:continue
        if len(seen)>=30 or used>=20_000_000:break
        seen.add(contig)
        try:
            ap=fetch('https://www.ebi.ac.uk/metagenomics/api/v2/assemblies/'+loc['assembly_accession']+'/analyses?page_size=100',folder/(loc['assembly_accession']+'.json'),1000000)
            data=json.loads(ap.read_text());analyses=data.get('items',data.get('results',[]))
            record=None
            for a in analyses[:3]:
                aid=a.get('accession') or a.get('mgya_accession')
                url='https://www.ebi.ac.uk/metagenomics/api/v1/analyses/'+aid+'/contigs/'+urllib.parse.quote(loc['contig_name'],safe='')
                try:
                    p=fetch(url,folder/(contig+'.fna'),min(5_000_000,20_000_000-used),headers={'Accept':'application/json'})
                    raw=p.read_bytes();used+=len(raw);record=SeqIO.read(io.StringIO(raw.decode()),'fasta');break
                except Exception as e:last_error=str(e)
            if record is None:raise ValueError(last_error if analyses else 'No analysis mapping')
            if record.id!=loc['contig_name']:raise ValueError('Returned contig differs')
            record.id=contig;record.annotations.update(molecule_type='DNA',topology='linear')
            dest=folder/(contig+'.gb');SeqIO.write(record,dest,'genbank')
            positions=locate_peptide(record,loc['peptide'],{'+':1,'-':-1}.get(loc['strand']))
            maps.append(dict(protein=loc['protein'],contig=contig,coordinates_json=json.dumps(positions),status='exact_unique_peptide_mapping' if len(positions)==1 else 'mapping_ambiguous_or_failed'))
            if len(positions)==1:sources[contig]=dict(accession=contig,path=dest.relative_to(ROOT).as_posix(),format='genbank',file_sha256=sha(dest),
                sequence_sha256=hashlib.sha256(str(record.seq).upper().encode()).hexdigest(),source_url=url,
                source_version_kind='MGnify_release_snapshot',database_release='2026_07',role='step2_MGnify_occurrence',length_nt=len(record),bytes=len(raw))
        except Exception as e:failures.append(dict(protein=loc['protein'],contig=contig,status='source_mapping_unresolved',reason=str(e)))
    table(folder/'acquired_source_records.tsv',list(sources.values()))
    table(folder/'peptide_mappings.tsv',maps);table(folder/'failures.tsv',failures)
    summary=dict(representatives_queried=len(reps),listed_occurrences=len(locations),previously_unacquired_contigs_attempted=len(seen),
        new_mapped_contigs=len(sources),new_DNA_bytes=used,failures=len(failures),
        limitation='Public API exposes representatives and listed occurrences only; complete cluster-member expansion requires separate bulk/Parquet lookup')
    write_json(folder/'summary.json',summary)
    return list(sources.values()),summary


from collections import Counter
