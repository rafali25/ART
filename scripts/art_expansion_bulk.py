"""Stream large DNA catalogs, deduplicate, call viral genes and export mapped pairs."""
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import asdict
import csv
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile
import time

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation, BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord
import pyrodigal_gv

from art_expansion_workflow import ROOT, OUT, OLD, CPUS, PROFILES, DATASETS, sha, file_md5, table, read_tsv, write_json, run
from art_reconstruction_common import parse_domains
from art_discovery import neighbors
from extract_art_loci import plan_locus,extract_record,write_gff


def dna_hash(dna):
    return min(hashlib.sha256(dna.encode()).hexdigest(),hashlib.sha256(str(Seq(dna).reverse_complement()).encode()).hexdigest())


def initialize_registry():
    conn=sqlite3.connect(OUT/'sequence_registry.sqlite')
    conn.execute('CREATE TABLE IF NOT EXISTS sequences (hash TEXT PRIMARY KEY, dataset TEXT, genome_key TEXT)')
    if conn.execute('SELECT count(*) FROM sequences').fetchone()[0]:return conn
    paths={}
    for folder in ['art_discovery','art_type1_public','art_type1_mgnify']:
        for r in read_tsv(ROOT/'results'/folder/'source_records.tsv'):
            if r.get('path') and r.get('format'):paths[r['path']]=r['format']
    for r in read_tsv(OLD/'panel_acquisition.tsv'):
        if r.get('status')=='retrieved':paths[r['path']]=r['format']
    for sub in ['deferred_acquisition/processed','diverse_partner_search']:
        for r in read_tsv(OLD/sub/'acquired_source_records.tsv'):paths[r['path']]=r['format']
    for name,fmt in paths.items():
        path=ROOT/name.replace('\\','/')
        if path.is_file():
            for r in SeqIO.parse(path,fmt):conn.execute('INSERT OR IGNORE INTO sequences VALUES(?,?,?)',(dna_hash(str(r.seq).upper()),'previous_workspace',r.id))
    jug=OLD/'jug_screen/sources/deposited/139.curated.Mac.phage.genomes.fa'
    if jug.exists():
        for r in SeqIO.parse(jug,'fasta'):conn.execute('INSERT OR IGNORE INTO sequences VALUES(?,?,?)',(dna_hash(str(r.seq).upper()),'Jug_pilot',r.id))
    conn.commit();return conn


def gene_job(item):
    key,original_id,dna=item
    finder=pyrodigal_gv.ViralGeneFinder(meta=True)
    called=finder.find_genes(dna.encode());rows=[];codes=Counter()
    for i,g in enumerate(called,1):
        a,b=g.begin-1,g.end
        cds=Seq(dna[a:b]);cds=cds if g.strand==1 else cds.reverse_complement()
        text=g.translate().rstrip('*');actual=str(cds.translate(table=g.translation_table)).rstrip('*')
        p5=bool(g.partial_begin if g.strand==1 else g.partial_end)
        p3=bool(g.partial_end if g.strand==1 else g.partial_begin)
        if not p5 and actual:actual='M'+actual[1:]
        compatible=len(text)==len(actual) and all(a==b or b=='X' for a,b in zip(actual,text))
        if not compatible:raise ValueError(f'Unexplained gene translation mismatch: {key} {i}')
        codes[g.translation_table]+=1
        rows.append((key+f'_pgv{i}',key,a,b,g.strand,g.translation_table,int(p5),int(p3),text,
            hashlib.sha256(str(cds).encode()).hexdigest(),hashlib.sha256(text.encode()).hexdigest(),int(actual==text)))
    return rows,str(called.metagenomic_bin),dict(codes)


def bulk_screen(name):
    folder=OUT/name;folder.mkdir(exist_ok=True)
    summarypath=folder/'summary.json'
    if summarypath.exists() and json.loads(summarypath.read_text()).get('screen_complete'):
        print(name+' already screened; retaining frozen output',flush=True);return json.loads(summarypath.read_text())
    started=time.monotonic();source=OUT/'downloads'/DATASETS[name]['file']
    assert source.exists()
    assert file_md5(source)==DATASETS[name]['md5'],'Catalog checksum changed'
    registry=initialize_registry()
    with tempfile.TemporaryDirectory(prefix=f'art-{os.environ["SLURM_JOB_ID"]}-{name}-',dir='/tmp') as temp:
        work=Path(temp);db=sqlite3.connect(work/'genes.sqlite');db.row_factory=sqlite3.Row
        db.execute('CREATE TABLE genomes (key TEXT PRIMARY KEY, original_id TEXT, header TEXT, dna TEXT, model TEXT, codes TEXT)')
        db.execute('CREATE TABLE genes (key TEXT PRIMARY KEY, genome_key TEXT, start0 INT,end0 INT,strand INT,code INT,partial5 INT,partial3 INT,seq TEXT,CDS_sha TEXT,protein_sha TEXT,exact INT)')
        fasta=work/'proteins.faa';genesgz=folder/'gene_coordinates.tsv.gz'
        genome_rows=[];pending=[];total=0;skipped=0;protein_count=0;exact_diffs=0;codecounts=Counter();models=Counter()
        with gzip.open(source,'rt') as inp, fasta.open('w') as aa, gzip.open(genesgz,'wt') as gf, ProcessPoolExecutor(max_workers=CPUS) as pool:
            gw=csv.writer(gf,delimiter='\t',lineterminator='\n')
            gw.writerow(['protein_key','genome_key','start0','end0','strand','genetic_code','partial_5prime','partial_3prime','protein_length_aa','CDS_sha256','protein_sha256','translation_exact'])
            def finish_batch(batch):
                nonlocal protein_count,exact_diffs
                jobs=[(x['key'],x['original_id'],x['dna']) for x in batch]
                for item,(gs,model,codes) in zip(batch,pool.map(gene_job,jobs)):
                    db.execute('INSERT INTO genomes VALUES(?,?,?,?,?,?)',(item['key'],item['original_id'],item['header'],item['dna'],model,json.dumps(codes)))
                    db.executemany('INSERT INTO genes VALUES(?,?,?,?,?,?,?,?,?,?,?,?)',gs)
                    for g in gs:
                        aa.write('>'+g[0]+' '+item['original_id']+'\n'+g[8]+'\n')
                        gw.writerow([*g[:8],len(g[8]),*g[9:]])
                    protein_count+=len(gs);exact_diffs+=sum(not g[-1] for g in gs)
                    codecounts.update({str(k):v for k,v in codes.items()});models[model]+=1
                    genome_rows.append(dict(genome_key=item['key'],original_id=item['original_id'],source_header=item['header'],length_nt=len(item['dna']),
                        canonical_sequence_sha256=item['hash'],status='called',genes=len(gs),selected_model=model,genetic_codes=json.dumps(codes)))
                db.commit()
                print(f'{name}: {total} records read, {protein_count} proteins, {skipped} duplicate/invalid records',flush=True)
            for index,r in enumerate(SeqIO.parse(inp,'fasta'),1):
                total+=1;dna=str(r.seq).upper();key=name+f'{index:06d}';dh=dna_hash(dna)
                existing=registry.execute('SELECT dataset,genome_key FROM sequences WHERE hash=?',(dh,)).fetchone()
                if existing:
                    skipped+=1;genome_rows.append(dict(genome_key=key,original_id=r.id,source_header=r.description,length_nt=len(dna),canonical_sequence_sha256=dh,
                        status='exact_DNA_duplicate',duplicate_dataset=existing[0],duplicate_genome_key=existing[1]));continue
                if not set(dna)<=set('ACGTRYSWKMBDHVN'):
                    skipped+=1;genome_rows.append(dict(genome_key=key,original_id=r.id,status='invalid_DNA_symbols',length_nt=len(dna)));continue
                registry.execute('INSERT INTO sequences VALUES(?,?,?)',(dh,name,key))
                pending.append(dict(key=key,original_id=r.id,header=r.description,dna=dna,hash=dh))
                if len(pending)==64:finish_batch(pending);pending=[]
            if pending:finish_batch(pending)
        db.execute('CREATE INDEX genomic_order ON genes(genome_key,start0,end0)');db.commit()
        table(folder/'genome_manifest.tsv',genome_rows)
        if not protein_count:
            for name_of_table in ['profile_hits.tsv','pairs.tsv','acquired_source_records.tsv','RT_review_hits.tsv']:
                table(folder/name_of_table,[])
            for name_of_fasta in ['RT_proteins.faa','partners.faa','RT_cores.faa']:
                (folder/name_of_fasta).write_text('')
            db.close();registry.commit();registry.close()
            summary=dict(dataset=name,records_in_archive=total,unique_new_genomes=0,duplicate_or_invalid_records=skipped,
                proteins_searched=0,primary_pairs=0,complete_passing_typeI_pairs=0,complete_passing_typeII_pairs=0,
                screen_complete=True,source_archive_sha256=sha(source),source_archive_md5=DATASETS[name]['md5'],
                interpretation='No new callable DNA remained after exact deduplication; no empty-database HMM search attempted')
            write_json(summarypath,summary);return summary
        # The gene database and large temporary protein FASTA stay on node-local
        # scratch. Coordinate/protein hashes plus original DNA and caller settings
        # retain reproducibility; only selected sequences need durable FASTAs.
        search=folder/'searches';search.mkdir(exist_ok=True)
        def search_one(item):
            family,profile=item;stem=search/family
            run([ROOT/'.venv/bin/hmmsearch','--cpu','1','--seed','42','--noali','-Z',protein_count,'--domZ',protein_count,
                '-E','.001','--domE','.001','--domtblout',str(stem)+'.domtbl',profile,fasta],str(stem)+'.log')
            return [dict(h,family=family) for h in parse_domains(str(stem)+'.domtbl')]
        hits=[]
        with ThreadPoolExecutor(max_workers=min(4,CPUS)) as executor:
            for hs in executor.map(search_one,PROFILES.items()):hits.extend(hs)
        for h in hits:h['passes_frozen_screen']=h['sequence_evalue']<=1e-5 and h['domain_evalue']<=1e-5 and h['profile_coverage']>=.6
        table(folder/'profile_hits.tsv',hits)
        result=export_pairs(db,hits,folder,name)
        db.close();registry.commit();registry.close()
        summary=dict(dataset=name,records_in_archive=total,unique_new_genomes=total-skipped,duplicate_or_invalid_records=skipped,
            proteins_searched=protein_count,reported_domains=len(hits),passing_domains=dict(Counter(h['family'] for h in hits if h['passes_frozen_screen'])),
            genetic_codes=dict(codecounts),caller_models=dict(models),ambiguous_codon_exact_differences=exact_diffs,
            unexplained_translation_errors=0,**result,seconds=round(time.monotonic()-started,3),
            source_archive_sha256=sha(source),source_archive_md5=DATASETS[name]['md5'],profiles={f:sha(p) for f,p in PROFILES.items()},
            caller_version=pyrodigal_gv.__version__,search_space_Z=protein_count,screen_complete=True,
            source_annotation='source DNA only; new viral gene calls; no source CDS annotations supplied',
            array_significance_tests_performed=False,interpretation='Mapped computational candidates; no functional ART validation')
        write_json(summarypath,summary)
        print(name+' DONE: '+json.dumps(summary),flush=True)
        return summary


def export_pairs(db,hits,folder,name):
    best=defaultdict(dict)
    for h in hits:
        if h['passes_frozen_screen'] or (h['family']=='RVT_1' and h['domain_evalue']<=.001):
            old=best[h['family']].get(h['target'])
            if old is None or (h['family']=='RVT_1' and h['protein_envelope_start0']<old['protein_envelope_start0']) or (h['family']!='RVT_1' and h['domain_bits']>old['domain_bits']):best[h['family']][h['target']]=h
    partners={key:family for family in ['TypeI_two_seed','Type_I_partner','Type_II_partner'] for key in best[family]}
    rtkeys=set(best['ART_RT'])|set(best['RVT_1'])
    candidates=[];gene_cache={};genome_cache={}
    def genome_genes(gkey):
        if gkey not in gene_cache:
            gene_cache[gkey]=[dict(r,gene_key=r['key'],protein_key=r['key'],accession=gkey) for r in db.execute('SELECT * FROM genes WHERE genome_key=? ORDER BY start0,end0',(gkey,))]
        return gene_cache[gkey]
    for pkey,family in partners.items():
        p=dict(db.execute('SELECT * FROM genes WHERE key=?',(pkey,)).fetchone());p['gene_key']=pkey
        genes=genome_genes(p['genome_key'])
        upstream=neighbors(p,genes,'upstream',3)
        for rank,g in enumerate(upstream,1):
            if g['key'] not in rtkeys or g['strand']!=p['strand']:continue
            if any(v['key'] in rtkeys for v in upstream[:rank-1]):continue
            downstream=neighbors(g,genes,'downstream',3)
            if pkey not in {v['key'] for v in downstream}:continue
            candidates.append(dict(rt=g,partner=p,family='Type_II' if family=='Type_II_partner' else 'Type_I',
                partner_profile=family,upstream_gene_rank=rank,partner_downstream_rank=next(i for i,v in enumerate(downstream,1) if v['key']==pkey)))
    byrt=defaultdict(list)
    for p in candidates:byrt[p['rt']['key']].append(p)
    chosen=[]
    for key,ps in byrt.items():
        ps.sort(key=lambda p:(p['partner_downstream_rank'],p['upstream_gene_rank'],p['partner']['key']))
        for i,p in enumerate(ps):p['association']='primary_nearest_compatible' if i==0 else 'alternative';chosen.append(p)
    pairrows=[];rtrecords=[];partnerrecords=[];corerecords=[];source_manifest=[];emitted=set();source_dir=folder/'sources';source_dir.mkdir(exist_ok=True)
    locus_dir=folder/'loci';locus_dir.mkdir(exist_ok=True)
    for p in chosen:
        g,pg=p['rt'],p['partner'];gkey=g['genome_key'];lid=gkey+'__'+g['key'];rtgate=g['key'] in best['ART_RT']
        if p['association']!='primary_nearest_compatible':lid+='__alternative_'+pg['key']
        original=dict(db.execute('SELECT * FROM genomes WHERE key=?',(gkey,)).fetchone())
        rec=SeqRecord(Seq(original['dna']),id=gkey,name=gkey,description=original['header'])
        rec.annotations.update(molecule_type='DNA',topology='linear',comment='Source topology not established; no boundary wrapping enabled.')
        features={}
        for v in genome_genes(gkey):
            a=BeforePosition(v['start0']) if (v['partial5'] if v['strand']==1 else v['partial3']) else v['start0']
            b=AfterPosition(v['end0']) if (v['partial3'] if v['strand']==1 else v['partial5']) else v['end0']
            f=SeqFeature(SimpleLocation(a,b,strand=v['strand']),type='CDS',qualifiers={'protein_id':[v['key']], 'translation':[v['seq']], 'transl_table':[str(v['code'])], 'inference':['ab initio prediction:Pyrodigal-gv:'+pyrodigal_gv.__version__]})
            rec.features.append(f);features[v['key']]=f
        src=source_dir/(gkey+'.gb')
        if gkey not in emitted:
            SeqIO.write(rec,src,'genbank');emitted.add(gkey)
            source_manifest.append(dict(accession=gkey,original_id=original['original_id'],source_header=original['header'],path=src.relative_to(ROOT).as_posix(),format='genbank',
                file_sha256=sha(src),sequence_sha256=hashlib.sha256(original['dna'].encode()).hexdigest(),role=name+'_public_catalog_candidate',
                source_version_kind='published_DNA_catalog_snapshot',database_release=name+'_pinned_archive_MD5_'+DATASETS[name]['md5']))
        plan=plan_locus(rec,features[g['key']],[features[pg['key']]])
        locus=extract_record(rec,features[g['key']],plan,lid)
        SeqIO.write(locus,locus_dir/(lid+'.gb'),'genbank');SeqIO.write(locus,locus_dir/(lid+'.fna'),'fasta');write_gff(locus_dir/(lid+'.gff3'),locus)
        write_json(locus_dir/(lid+'.map.json'),dict(source_genome_key=gkey,original_id=original['original_id'],source_path=src.relative_to(ROOT).as_posix(),
            segments=[asdict(x) for x in plan['segments']],upstream_available_nt=plan['upstream_available_nt']))
        rf=features[g['key']];pf=features[pg['key']]
        complete_checks=[]
        for f,v in [(rf,g),(pf,pg)]:
            cds=f.extract(rec.seq);computed=str(cds.translate(table=v['code'])).rstrip('*')
            if not v['partial5'] and computed:computed='M'+computed[1:]
            assert len(computed)==len(v['seq']) and all(a==b or b=='X' for a,b in zip(computed,v['seq']))
            try:complete_checks.append(not(v['partial5'] or v['partial3'] or 'X' in v['seq']) and str(cds.translate(table=v['code'],cds=True))==v['seq'])
            except Exception:complete_checks.append(False)
        complete=all(complete_checks)
        rh=best['ART_RT'].get(g['key']);ph=best[p['partner_profile']][pg['key']];core=best['RVT_1'].get(g['key'])
        row=dict(locus_id=lid,dataset=name,original_source_id=original['original_id'],source_genome=gkey,RT_key=g['key'],partner_key=pg['key'],
            family=p['family'],association=p['association'],pair_complete=complete,ART_profile_passes=rtgate,
            RT_profile_coverage=rh['profile_coverage'] if rh else '',RT_profile_evalue=rh['sequence_evalue'] if rh else '',partner_profile=p['partner_profile'],
            partner_profile_coverage=ph['profile_coverage'],partner_evalue=ph['sequence_evalue'],
            RT_start0=g['start0'],RT_end0=g['end0'],partner_start0=pg['start0'],partner_end0=pg['end0'],strand=g['strand'],
            RT_length_aa=len(g['seq']),partner_length_aa=len(pg['seq']),RT_code=g['code'],partner_code=pg['code'],upstream_nt=plan['upstream_available_nt'],
            core_start0=max(0,core['protein_envelope_start0']-40) if core else '',core_end0=min(len(g['seq']),core['protein_envelope_end0']+40) if core else '',
            status='passing_profile_context_candidate' if rtgate else 'generic_RT_partner_context_requires_ART_family_review',array_status='not_assessed_in_new_catalog',
            source_path=src.relative_to(ROOT).as_posix(),locus_genbank=(locus_dir/(lid+'.gb')).relative_to(ROOT).as_posix(),coordinate_map=(locus_dir/(lid+'.map.json')).relative_to(ROOT).as_posix())
        pairrows.append(row)
        if p['association']=='primary_nearest_compatible':
            rtrecords.append(SeqRecord(Seq(g['seq']),id=lid,description=g['key']));partnerrecords.append(SeqRecord(Seq(pg['seq']),id=lid,description=pg['key']))
            if core:corerecords.append(SeqRecord(Seq(g['seq'][row['core_start0']:row['core_end0']]),id=lid,description='RVT_1_envelope_plus_minus_40aa'))
    table(folder/'pairs.tsv',pairrows);table(folder/'acquired_source_records.tsv',source_manifest)
    SeqIO.write(rtrecords,folder/'RT_proteins.faa','fasta');SeqIO.write(partnerrecords,folder/'partners.faa','fasta');SeqIO.write(corerecords,folder/'RT_cores.faa','fasta')
    # Retain every reported RT-family hit, including coverage failures, as a
    # separate review table; a partner absence does not make a NART reference.
    rtrows=[]
    for h in hits:
        if h['family'] in ['ART_RT','RVT_1']:
            v=dict(db.execute('SELECT * FROM genes WHERE key=?',(h['target'],)).fetchone())
            rtrows.append(dict(h,genome_key=v['genome_key'],start0=v['start0'],end0=v['end0'],strand=v['strand'],length_aa=len(v['seq'])))
    table(folder/'RT_review_hits.tsv',rtrows)
    prim=[r for r in pairrows if r['association']=='primary_nearest_compatible']
    return dict(mapped_pair_associations=len(pairrows),primary_pairs=len(prim),
        complete_passing_typeI_pairs=sum(r['family']=='Type_I' and r['pair_complete'] and r['ART_profile_passes'] for r in prim),
        complete_passing_typeII_pairs=sum(r['family']=='Type_II' and r['pair_complete'] and r['ART_profile_passes'] for r in prim),
        generic_RT_review_pairs=sum(not r['ART_profile_passes'] for r in prim),
        exported_candidate_genomes=len(emitted),pair_DNA_translation_checks_passed=True)
