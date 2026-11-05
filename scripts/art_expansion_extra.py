"""A bounded, source-linked UHGV stratum after the two jumbo-phage collections."""
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import io
import json
import math
import re

from Bio import SeqIO
import zstandard

from art_expansion_workflow import ROOT,OUT,DATASETS,fetch,file_md5,sha,table,write_json


def select_uhgv():
    from art_expansion_bulk import bulk_screen
    folder=OUT/'UHGV_selected';folder.mkdir(exist_ok=True)
    meta=fetch('https://zenodo.org/records/17402089/files/votus_metadata_extended.tsv.zst?download=1',
        folder/'votus_metadata_extended.tsv.zst',25_000_000,expected_md5='6702355349e933fb605e62a8bbaa50a9')
    decompressor=zstandard.ZstdDecompressor().stream_reader(meta.open('rb'))
    candidates=[]
    with io.TextIOWrapper(decompressor) as f:
        reader=csv.DictReader(f,delimiter='\t');schema=reader.fieldnames
        for r in reader:
            ids=[v for v in r.values() if v and re.fullmatch(r'UHGV-\d+',v)]
            if not ids:continue
            lengths=[(k,v) for k,v in r.items() if re.search(r'(genome|sequence).*(length|size)|^(length|size|length_nt)$',k,re.I) and v and re.fullmatch(r'\d+(?:\.0)?',v)]
            if not lengths:continue
            length=int(float(lengths[0][1]))
            if not 150_000<=length<=1_000_000:continue
            quality=';'.join(v for k,v in r.items() if re.search('quality|completeness',k,re.I) and v)
            taxonomy=';'.join(v for k,v in r.items() if re.search('tax|host|genus|family',k,re.I) and v)
            group=';'.join(v for k,v in r.items() if re.search('genus|family',k,re.I) and v) or 'unclassified'
            candidates.append(dict(genome_id=ids[0],length_nt=length,quality=quality,taxonomy=taxonomy,group=group,metadata_json=json.dumps(r)))
    write_json(folder/'metadata_schema.json',dict(columns=schema,long_genome_candidates=len(candidates),
        rule='150kb_to_1Mb; quality and Firmicutes-host annotations prioritize; source strata round-robin; maximum120 genomes and60MB planned transfers'))
    candidates.sort(key=lambda r:(not bool(re.search('Firmicutes|Bacill|Listeria|Staphyl',r['taxonomy'],re.I)),
        not bool(re.search('high|complete|100',r['quality'],re.I)),-r['length_nt'],r['genome_id']))
    groups=defaultdict(list)
    for r in candidates:groups[r['group']].append(r)
    chosen=[];reserved=0
    while len(chosen)<120 and any(groups.values()):
        for group in list(groups):
            if not groups[group] or len(chosen)>=120:continue
            r=groups[group].pop(0);cap=math.ceil(r['length_nt']*1.1)+10000
            if reserved+cap>60_000_000:continue
            r['transfer_cap']=cap;reserved+=cap;chosen.append(r)
    table(folder/'frozen_selection.tsv',chosen)
    if not chosen:
        raise ValueError('UHGV metadata schema produced no selectable records; inspect metadata_schema.json')
    source=OUT/'downloads/UHGV_selected.fasta.gz'
    if not source.exists():
        def acquire(r):
            gid=r['genome_id'];url=f'https://portal.nersc.gov/UHGV/votu_reps/{gid[:8]}/{gid}/{gid}.fna'
            try:
                path=fetch(url,folder/'sources'/(gid+'.fna'),r['transfer_cap'])
                rec=SeqIO.read(path,'fasta')
                if len(rec)!=r['length_nt']:raise ValueError('Genomic length differs from frozen metadata')
                return dict(r,status='retrieved',path=str(path.relative_to(ROOT)),source_url=url,sha256=sha(path)),rec
            except Exception as e:return dict(r,status='retrieval_failed',source_url=url,reason=str(e)),None
        outcomes=[];records=[]
        with ThreadPoolExecutor(max_workers=3) as executor:
            for result,rec in executor.map(acquire,chosen):
                outcomes.append(result)
                if rec:records.append(rec)
        table(folder/'retrieval_outcomes.tsv',outcomes)
        if not records:raise ValueError('No bounded UHGV source retrieval succeeded')
        with gzip.open(source,'wt') as f:SeqIO.write(records,f,'fasta')
    DATASETS['UHGV_selected']=dict(file=source.name,md5=file_md5(source))
    result=bulk_screen('UHGV_selected')
    write_json(folder/'selection_summary.json',dict(planned_genomes=len(chosen),planned_transfer_cap=reserved,
        public_metadata_source=str(meta.relative_to(ROOT)),screen=result,
        IMG_VR='Not additionally downloaded: UHGV already incorporates IMG/VR genomes; a separate authenticated/full IMG/VR expansion remains a future distinct cohort'))
    return json.loads((folder/'selection_summary.json').read_text())
