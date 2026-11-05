"""Bounded public mapped-genome pilot; architecture never defines sequence placement.

Internal and TSV coordinates are zero based, half open. Source annotations are
retained. Later stages are gated on validation plus measured cost and yield.
"""
from __future__ import annotations
import argparse
from collections import defaultdict,Counter
from dataclasses import asdict
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import time
import urllib.request

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition,BeforePosition,AfterPosition,SeqFeature,SimpleLocation
from Bio.SeqRecord import SeqRecord

from art_reconstruction_common import ROOT,Runner,sha,table,passes,stage_gate,validation_gate
from art_discovery import neighbors
from extract_art_loci import read_tsv,write_json,plan_locus,extract_record,pid,source_pieces,write_gff
from art_array_detector import refinement_scan,initial_scan,distinct_arrays,BASES
from art_array_assessment import aligned_pwm

FIELDS=['locus_id','source_accession_and_version','contig_id','rt_coordinates_and_strand',
 'art_profile_score_and_coverage','rt_completeness','ntd_boundary','partner_family_and_coordinates',
 'array_assessment_status','repeat_coordinates_and_consensus','spacer_lengths','upstream_sequence_available',
 'phylogenetic_assignment','discovery_route','final_evidence_category','protein_key','RT_protein_id',
 'record_role','source_path','source_file_sha256','oriented_locus_path','coordinate_map_path','notes']


def fetch_record(spec,out):
    acc=spec['accession'];fmt='embl' if spec['source']=='ENA' else 'genbank'
    url=('https://www.ebi.ac.uk/ena/browser/api/embl/'+acc+'?download=false') if fmt=='embl' else (
        'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+acc+'&rettype=gbwithparts&retmode=text')
    path=out/'sources'/(acc+('.embl' if fmt=='embl' else '.gb'));cached=path.exists()
    if not path.exists():
        with urllib.request.urlopen(url,timeout=55) as response:data=response.read(20000001)
        if len(data)>20000000:raise ValueError('Pilot per-record download limit exceeded')
        path.write_bytes(data)
    record=SeqIO.read(path,fmt)
    if record.id!=acc:raise ValueError(f'Accession/version mismatch: requested {acc}; returned {record.id}')
    if not any(f.type=='CDS' for f in record.features):raise ValueError('Stage 1 requires annotated CDSs: '+acc)
    return record,dict(accession=acc,path=path.relative_to(ROOT).as_posix(),format=fmt,record_role='public_stage1_pilot',
        selection_role=spec['selection_role'],source_url=url,length_nt=len(record),file_sha256=sha(path),
        sequence_sha256=hashlib.sha256(str(record.seq).upper().encode()).hexdigest(),retrieved_utc=datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat(),cache_reused=cached)


def load_assembly(row,out):
    """Validate real nucleotide/annotation linkage, optionally adding viral calls."""
    fasta=ROOT/row['nucleotide_path'];source=ROOT/row['annotation_path']
    if not fasta.is_file() or not source.is_file():raise ValueError('Missing nucleotide or annotation source for '+row['contig_id'])
    rec=SeqIO.read(source,row['annotation_format']);nt=SeqIO.read(fasta,'fasta')
    if rec.id!=row['source_accession_and_version'] or rec.id!=row['contig_id'] or nt.id!=rec.id:
        raise ValueError('Assembly accession/version or contig mapping mismatch')
    if str(nt.seq).upper()!=str(rec.seq).upper():raise ValueError('Annotation and nucleotide sequences differ')
    path=source;method=row['annotation_method'];gene_call_details={}
    if method=='pyrodigal_gv':
        if any(f.type=='CDS' for f in rec.features):raise ValueError('Do not overwrite existing CDS annotations with viral gene calls')
        import pyrodigal_gv,pyrodigal
        predicted=pyrodigal_gv.ViralGeneFinder(meta=True).find_genes(bytes(rec.seq))
        original_feature_count=len(rec.features)
        for i,gene in enumerate(predicted,1):
            # Prodigal: 1-based inclusive; Biopython: 0-based half-open.
            a=BeforePosition(gene.begin-1) if gene.partial_begin else ExactPosition(gene.begin-1)
            b=AfterPosition(gene.end) if gene.partial_end else ExactPosition(gene.end)
            rec.features.append(SeqFeature(SimpleLocation(a,b,strand=gene.strand),type='CDS',qualifiers={
                'protein_id':[rec.id+f'_pgv{i:05d}'],'locus_tag':[rec.id+f'_pgv{i:05d}'],
                'translation':[gene.translate().rstrip('*')],'transl_table':[str(gene.translation_table)],'codon_start':['1'],
                'product':['predicted viral protein'],'inference':['ab initio prediction:pyrodigal-gv:'+pyrodigal_gv.__version__],
                'note':['Original record retained; partial CDS requires review' if gene.partial_begin or gene.partial_end else 'Viral-model candidate gene call']}))
        folder=out/'derived_gene_calls';folder.mkdir(exist_ok=True);path=folder/(rec.id+'.gb')
        SeqIO.write(rec,path,'genbank')
        gff_record=rec[:]
        import copy
        gff_record.features=copy.deepcopy(rec.features)
        for index,feature in enumerate(gff_record.features):feature.qualifiers['source_feature_index']=[str(index)]
        write_gff(folder/(rec.id+'.gff3'),gff_record)
        gff_path=folder/(rec.id+'.gff3');lines=[]
        for line in gff_path.read_text().splitlines():
            fields=line.split('\t')
            if len(fields)==9 and 'pyrodigal-gv' in fields[8]:fields[1]='pyrodigal-gv';line='\t'.join(fields)
            lines.append(line)
        gff_path.write_text('\n'.join(lines)+'\n')
        gene_call_details=dict(pyrodigal_gv_version=pyrodigal_gv.__version__,pyrodigal_version=pyrodigal.__version__,
            mode='ViralGeneFinder(meta=True); standard plus giant-virus/alternative-code models',
            selected_model=str(predicted.metagenomic_bin),original_features_preserved=original_feature_count,
            predicted_CDS=len(predicted),genetic_codes=dict(Counter(g.translation_table for g in predicted)))
    elif method!='source_original':raise ValueError('Unrecognized annotation method: '+method)
    if not any(f.type=='CDS' for f in rec.features):raise ValueError('No annotated or predicted genes in assembly')
    meta=dict(accession=rec.id,path=path.relative_to(ROOT).as_posix(),format='genbank' if method=='pyrodigal_gv' else row['annotation_format'],
        record_role='public_stage2_pilot',source_url=row['source_url'],selection_role=row['selection_reason'],length_nt=len(rec),
        file_sha256=sha(path),sequence_sha256=hashlib.sha256(str(rec.seq).upper().encode()).hexdigest(),
        annotation_method=method,original_annotation_path=row['annotation_path'],original_annotation_sha256=sha(source),
        nucleotide_path=row['nucleotide_path'],nucleotide_file_sha256=sha(fasta),gene_call_details=json.dumps(gene_call_details),
        source_record_status='UNVERIFIED' if 'UNVERIFIED' in rec.description else 'deposited_sequence')
    return rec,meta


def ingest(records,out):
    genes={};proteins=[];features={};by_contig=defaultdict(list);issues=[]
    for rec,meta in records:
        if '.' not in rec.id or not str(rec.seq):raise ValueError('Versioned retrievable nucleotide record required')
        for i,f in enumerate(rec.features):
            if f.type!='CDS' or not f.location:continue
            translation=f.qualifiers.get('translation',[''])[0].rstrip('*')
            if not translation:
                issues.append(dict(contig_id=rec.id,feature_index=i,status='CDS_without_translation_not_searched'));continue
            if set(translation)-set('ACDEFGHIKLMNPQRSTVWYXBZJUO'):
                issues.append(dict(contig_id=rec.id,feature_index=i,status='ambiguous_or_internal_stop_protein_not_searched'));continue
            ambiguous_residues=sum(c not in 'ACDEFGHIKLMNPQRSTVWY' for c in translation)
            if ambiguous_residues:issues.append(dict(contig_id=rec.id,feature_index=i,status='ambiguous_amino_acids_retained_for_HMM_search_and_flagged'))
            key='g'+hashlib.sha256((rec.id+'|'+str(i)+'|'+str(f.location)).encode()).hexdigest()[:16]
            coding=f.extract(rec.seq);table_id=int(f.qualifiers.get('transl_table',['11'])[0]);codon_start=int(f.qualifiers.get('codon_start',['1'])[0])
            fuzzy=any(not isinstance(p.start,ExactPosition) or not isinstance(p.end,ExactPosition) for p in f.location.parts)
            try:
                computed=str(coding[codon_start-1:].translate(table=table_id,cds=True));transqc='complete_CDS_translation_matches' if computed==translation else 'translation_disagrees'
            except Exception:transqc='complete_CDS_translation_not_confirmed'
            g=dict(protein_key=key,gene_key=key,accession=rec.id,contig_id=rec.id,protein_id=pid(f) or f'CDS_{i}',
                start0=int(f.location.start),end0=int(f.location.end),strand=f.location.strand,feature_index=i,
                source_parts_json=json.dumps([[int(p.start),int(p.end),p.strand] for p in f.location.parts]),
                product=f.qualifiers.get('product',[''])[0],length_aa=len(translation),contig_length_nt=len(rec),
                translation_qc=transqc,fuzzy_boundary=fuzzy,codon_start=codon_start,genetic_code=table_id,
                ambiguous_residue_count=ambiguous_residues,
                annotation_method=meta.get('annotation_method','source_original'),source_record_status=meta.get('source_record_status','deposited_sequence'),
                source_path=meta['path'],source_sha256=meta['file_sha256'],record_role=meta['record_role'])
            genes[key]=g;features[key]=f;by_contig[rec.id].append(g)
            proteins.append(SeqRecord(Seq(translation),id=key,description=f'{rec.id} {g["protein_id"]}'))
    SeqIO.write(proteins,out/'proteins.faa','fasta');table(out/'protein_to_contig.tsv',list(genes.values()))
    table(out/'annotation_issues.tsv',issues,fields=['contig_id','feature_index','status'])
    return genes,features,by_contig,proteins


def find_candidates(genes,by_contig,hits,cfg):
    art={};rvt={};partner={};all_art={}
    for h in hits:
        key=h['target']
        if h['search_family']=='ART_RT':
            if key not in all_art or h['domain_bits']>all_art[key]['domain_bits']:all_art[key]=h
            if passes(h,cfg) and (key not in art or h['domain_bits']>art[key]['domain_bits']):art[key]=h
        elif h['search_family']=='RVT_1' and h['domain_evalue']<=.001:
            if key not in rvt or h['protein_envelope_start0']<rvt[key]['protein_envelope_start0']:rvt[key]=h
        elif h['search_family']!='RVT_1' and passes(h,cfg):partner.setdefault(key,[]).append(h)
    routes=defaultdict(set)
    for key in art:routes[key].add('ART_protein_first')
    orphan=[];upstream_rows=[]
    for key,hs in partner.items():
        g=genes[key];up=neighbors(g,by_contig[g['accession']],'upstream',3);found=False
        for rank,other in enumerate(up,1):
            upstream_rows.append(dict(partner_key=key,partner_family=';'.join(h['search_family'] for h in hs),
                upstream_rank=rank,**{k:other[k] for k in ('protein_key','accession','protein_id','start0','end0','strand')},
                same_orientation=other['strand']==g['strand'],ART_profile_match=other['protein_key'] in all_art,
                generic_RT_match=other['protein_key'] in rvt))
            if other['protein_key'] in rvt or other['protein_key'] in all_art:
                routes[other['protein_key']].add('partner_first');found=True
        if not found:orphan.append(dict(partner_key=key,accession=g['accession'],partner_families=';'.join(h['search_family'] for h in hs),status='partner_homology_without_RT_in_first_three_upstream_genes'))
    return routes,art,rvt,partner,all_art,upstream_rows,orphan


def array_architecture(record,rt,plan,locus,out,cfg,baseline=None):
    rts=next(f for f in record.features if f.type=='CDS' and pid(f)==pid(rt));start=int(rts.location.start)
    upstream=str(record.seq[max(0,start-6000):start]).upper();offset=max(0,start-6000)
    if baseline:
        rr=[r for r in read_tsv(ROOT/'results/art_arrays_assessed/repeats.tsv') if r['array_id']==baseline['array_id']]
        repeats=[dict(locus_start0=int(r['locus_start0']),locus_end0=int(r['locus_end0']),sequence=r['sequence'],
            source_segments=json.loads(r['source_segments_json'])) for r in rr]
        return baseline['assessment_status'],repeats,json.loads(baseline['spacer_lengths_nt_json']),dict(provenance='previous_nine_locus_assessment',array_id=baseline['array_id'])
    scan=refinement_scan(upstream,cfg);initial=initial_scan(upstream[-3000:],cfg)
    candidates=distinct_arrays(scan['candidates']);model={};repeats=[];spacers=[]
    if candidates:
        model=aligned_pwm(upstream,candidates[0],cfg,out,locus)
        intervals=model.get('final_intervals',[])
        for a,b in intervals:
            repeats.append(dict(locus_start0=a+offset,locus_end0=b+offset,sequence=upstream[a:b],
                source_segments=source_pieces(a+offset,b+offset,plan['segments'])))
        spacers=[c-b for (a,b),(c,d) in zip(intervals,intervals[1:])]
    evidence=dict(provenance='new_exploratory_scan_not_calibrated_for_acceptance',initial=initial,extended=scan,model=model,
        no_empirical_array_pvalue_claim=True)
    write_json(out/'array_details'/(locus+'.json'),evidence)
    return 'exploratory_array_requires_validation_and_cohort_nulls' if len(repeats)>=3 else 'no_supported_array_call_exploratory_scan',repeats,spacers,evidence


def run_stage1(cfg,out,stage=1):
    started=time.perf_counter();run=Runner(out);records=[];source_rows=[];seen=set();aliases=[]
    if cfg['include_baseline_for_reconstruction']:
        previous=out.parent/f'stage{stage-1}'/'source_records.tsv'
        sources=read_tsv(previous) if stage>1 else read_tsv(ROOT/'results/art_discovery/source_records.tsv')
        for r in sources:
            path=ROOT/r['path']
            if sha(path)!=r['file_sha256']:raise ValueError('Baseline source changed: '+r['path'])
            rec=SeqIO.read(path,r['format']);m=dict(r,record_role='baseline_reference' if stage==1 else r['record_role'],source_url=r.get('source_url','previously_archived_source'))
            records.append((rec,m));source_rows.append(m);seen.add(r['sequence_sha256'])
    specs=cfg['stage1_records'] if stage==1 else read_tsv(ROOT/cfg[f'stage{stage}_assembly_manifest'])
    limit=cfg['maximum_pilot_records'] if stage==1 else cfg[f'stage{stage}_maximum_records']
    if not specs or len(specs)>limit:raise ValueError('Empty panel or pilot record budget exceeded')
    for spec in specs:
        rec,m=fetch_record(spec,out) if stage==1 else load_assembly(spec,out)
        if stage>1:m['record_role']=f'public_stage{stage}_pilot'
        if m['sequence_sha256'] in seen:
            aliases.append(dict(accession=m['accession'],status='exact_sequence_duplicate_not_searched_twice'));continue
        records.append((rec,m));source_rows.append(m);seen.add(m['sequence_sha256'])
    new_role=f'public_stage{stage}_pilot'
    new_sources=[m for r,m in records if m['record_role']==new_role]
    if sum((ROOT/m['path']).stat().st_size for m in new_sources)>cfg['maximum_pilot_download_bytes']:raise ValueError('Pilot download budget exceeded')
    table(out/'source_records.tsv',source_rows);table(out/'source_aliases.tsv',aliases,fields=['accession','status'])
    genes,features,by_contig,proteins=ingest(records,out)
    all_hits=[];search_seconds=0.
    for family in ('ART_RT','Type_I_partner','Type_II_partner','RVT_1'):
        hmm=run.binaries['RVT_1'] if family=='RVT_1' else ROOT/'results/art_discovery/profiles'/(family+'.hmm')
        hs,elapsed=run.hmm(hmm,out/'proteins.faa',out/'searches'/family,z=len(proteins))
        for h in hs:h['search_family']=family
        all_hits+=hs;search_seconds+=elapsed
    table(out/'profile_hits.tsv',all_hits)
    routes,art,rvt,partners,all_art,up,orphans=find_candidates(genes,by_contig,all_hits,cfg)
    table(out/'partner_first_upstream_genes.tsv',up);table(out/'partner_orphans.tsv',orphans,fields=['partner_key','accession','partner_families','status'])
    record_map={r.id:r for r,m in records};meta_map={m['accession']:m for r,m in records}
    for orphan in orphans:
        g=genes[orphan['partner_key']];rec=record_map[g['accession']];a=max(0,g['start0']-6000);b=min(len(rec),g['end0']+6000)
        segment=rec[a:b];segment.id=g['protein_key'];SeqIO.write(segment,out/'loci'/(g['protein_key']+'.partner_window.gb'),'genbank')
    baseline={r['accession']+'__'+r['RT_protein_id']:r for r in read_tsv(ROOT/'results/art_arrays_assessed/array_assessments.tsv')}
    acfg=json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
    rows=[];architecture=[];seqrows=[];candidates=[]
    gate=validation_gate()
    for key in sorted(routes):
        g=genes[key];rec=record_map[g['accession']];rt=features[key];loc=g['accession']+'__'+g['protein_id']
        downstream=neighbors(g,by_contig[g['accession']],'downstream',3);ps=[];pfeatures=[]
        for rank,pg in enumerate(downstream,1):
            for hit in partners.get(pg['protein_key'],[]):
                ps.append(dict(family=hit['search_family'],protein_id=pg['protein_id'],start0=pg['start0'],end0=pg['end0'],
                    strand=pg['strand'],downstream_gene_rank=rank,compatible_adjacent_context=rank==1 and pg['strand']==g['strand'],
                    domain_bits=hit['domain_bits'],profile_coverage=hit['profile_coverage']))
                pfeatures.append(features[pg['protein_key']])
        plan=plan_locus(rec,rt,pfeatures);oriented=extract_record(rec,rt,plan,loc)
        SeqIO.write(oriented,out/'loci'/(loc+'.gb'),'genbank');SeqIO.write(oriented,out/'loci'/(loc+'.fna'),'fasta');write_gff(out/'loci'/(loc+'.gff3'),oriented)
        write_json(out/'loci'/(loc+'.map.json'),dict(segments=[asdict(s) for s in plan['segments']],
            metadata={k:v for k,v in plan.items() if k not in ('segments','plus_segments')},source_accession=rec.id,
            source_path=g['source_path'],source_sha256=g['source_sha256']))
        ast,reps,spacers,ae=array_architecture(oriented,rt,plan,loc,out,acfg,baseline.get(loc))
        consensus=''
        if reps and len({len(r['sequence']) for r in reps})==1:
            consensus=''.join(max(BASES,key=lambda b:sum(r['sequence'][j]==b for r in reps)) for j in range(len(reps[0]['sequence'])))
        rvt_hit=rvt.get(key);ntd=dict(status='RVT_1_not_detected') if not rvt_hit else dict(start0=0,end0=rvt_hit['protein_envelope_start0'],
            coordinate_system='protein_0_based_half_open',definition='residues_before_first_significant_PF00078_envelope',
            profile_accession=run.rvt_accession,domain_evalue=rvt_hit['domain_evalue'])
        ah=all_art.get(key);profile={} if not ah else dict(bits=ah['domain_bits'],coverage=ah['profile_coverage'],sequence_evalue=ah['sequence_evalue'],domain_evalue=ah['domain_evalue'])
        complete=dict(translation_qc=g['translation_qc'],fuzzy_boundary=g['fuzzy_boundary'],length_aa=g['length_aa'],
            ambiguous_residue_count=g['ambiguous_residue_count'],published_350_900_filter=350<=g['length_aa']<=900,
            gene_call_review=g['translation_qc']!='complete_CDS_translation_matches' or g['fuzzy_boundary'] or g['ambiguous_residue_count']>0)
        upstream_genes=neighbors(g,by_contig[g['accession']],'upstream',3)
        noncoding_gap=(g['start0']-upstream_genes[0]['end0'] if g['strand']==1 else upstream_genes[0]['start0']-g['end0']) if upstream_genes else plan['upstream_available_nt']
        row=dict(locus_id=loc,source_accession_and_version=rec.id,contig_id=rec.id,rt_coordinates_and_strand=json.dumps(dict(start0=g['start0'],end0=g['end0'],strand=g['strand'],parts=json.loads(g['source_parts_json']))),
            art_profile_score_and_coverage=json.dumps(profile),rt_completeness=json.dumps(complete),ntd_boundary=json.dumps(ntd),
            partner_family_and_coordinates=json.dumps(ps),array_assessment_status=ast,repeat_coordinates_and_consensus=json.dumps(dict(consensus=consensus,repeats=reps)),
            spacer_lengths=json.dumps(spacers),upstream_sequence_available=plan['upstream_available_nt'],phylogenetic_assignment='pending_sequence_only_domain_tree',
            discovery_route=';'.join(sorted(routes[key])),final_evidence_category='mapped_RT_candidate_requires_validation',protein_key=key,
            RT_protein_id=g['protein_id'],record_role=g['record_role'],source_path=g['source_path'],source_file_sha256=g['source_sha256'],
            oriented_locus_path=(out/'loci'/(loc+'.gb')).relative_to(ROOT).as_posix(),coordinate_map_path=(out/'loci'/(loc+'.map.json')).relative_to(ROOT).as_posix(),
            notes='No definitive ART/NART assignment; no complete-locus promotion while validation gate is incomplete; annotation='+g['annotation_method']+'; source_status='+g['source_record_status'])
        rows.append(row)
        architecture.append(dict(locus_id=loc,protein_key=key,partner_family_and_coordinates=row['partner_family_and_coordinates'],
            upstream_noncoding_gap_nt=max(0,noncoding_gap),array_assessment_status=ast,repeat_coordinates_and_consensus=row['repeat_coordinates_and_consensus'],
            spacer_lengths=row['spacer_lengths'],upstream_sequence_available=plan['upstream_available_nt'],
            architecture_typeI_noncoding500_anchor=any(p['family']=='Type_I_partner' and p['compatible_adjacent_context'] for p in ps) and noncoding_gap>=500))
        seqrows.append(dict(protein_key=key,locus_id=loc,source_accession_and_version=rec.id,contig_length_nt=len(rec),
            RT_protein_id=g['protein_id'],length_aa=g['length_aa'],ART_profile_score=ah['domain_bits'] if ah else '',
            ART_profile_coverage=ah['profile_coverage'] if ah else '',ART_screen_pass=key in art,
            RVT_1_envelope_start0=rvt_hit['protein_envelope_start0'] if rvt_hit else '',RVT_1_envelope_end0=rvt_hit['protein_envelope_end0'] if rvt_hit else '',
            ntd_boundary=json.dumps(ntd),placement='pending_sequence_only_domain_tree',record_role=g['record_role']))
        candidates.append(next(p for p in proteins if p.id==key))
    table(out/'candidate_loci.tsv',rows,FIELDS);table(out/'observed_locus_architecture.tsv',architecture);table(out/'sequence_family_placement.tsv',seqrows)
    SeqIO.write(candidates,out/'candidate_RT_proteins.faa','fasta')
    elapsed=time.perf_counter()-started
    metrics=dict(stage=stage,records_searched=len(new_sources),baseline_records_included=len(records)-len(new_sources),
        nucleotide_nt=sum(len(r) for r,m in records if m['record_role']==new_role),proteins_searched=len(proteins),
        new_proteins=sum(g['record_role']==new_role for g in genes.values()),search_seconds=round(search_seconds,3),elapsed_seconds=round(elapsed,3),
        new_RT_candidates=sum(r['record_role']==new_role for r in rows),total_RT_candidates=len(rows),
        new_partner_orphans=sum(genes[r['partner_key']]['record_role']==new_role for r in orphans),
        input_file_bytes=sum((ROOT/m['path']).stat().st_size for m in new_sources),new_candidate_protein_aa=sum(g['length_aa'] for k,g in genes.items() if k in routes and g['record_role']==new_role),
        selection_bias='purposive near-relative and comparison panel; yield is not a prevalence or sensitivity estimate',
        HMM_effective_database_size=len(proteins),available_partner_profiles=['Type_I_partner','Type_II_partner'],
        unavailable_partner_profiles=['Type_III_first_partner','Type_III_second_partner'],validation_gate=gate,
        next_stage_gate=stage_gate(stage+1,cfg,dict(records_searched=len(new_sources),search_seconds=search_seconds,new_RT_candidates=sum(r['record_role']==new_role for r in rows))))
    write_json(out/'cost_and_yield.json',metrics);write_json(out/'run_config.json',cfg)
    write_json(out/'source_provenance.json',dict(software=run.metadata,source_records=source_rows,config_sha256=sha(ROOT/'config/art_expansion.json')))
    print(json.dumps(metrics,indent=2))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',type=int,choices=(1,2,3),default=1);ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--config',type=Path,default=ROOT/'config/art_expansion.json');args=ap.parse_args()
    cfg=json.loads(args.config.read_text());base=ROOT/cfg['output'];out=base/f'stage{args.stage}'
    for d in ('sources','searches','loci','alignments','array_details'):(out/d).mkdir(parents=True,exist_ok=True)
    if args.stage!=1:
        prev=base/f'stage{args.stage-1}'/'cost_and_yield.json';metrics=json.loads(prev.read_text()) if prev.exists() else None
        gate=stage_gate(args.stage,cfg,metrics);write_json(out/'stage_gate.json',gate)
        if args.prepare_only:
            rows=read_tsv(ROOT/cfg[f'stage{args.stage}_assembly_manifest'])
            if len(rows)>cfg[f'stage{args.stage}_maximum_records']:raise ValueError('Assembly preparation budget exceeded')
            prepared=[load_assembly(r,out)[1] for r in rows]
            table(out/'prepared_assembly_sources.tsv',prepared)
            write_json(out/'preparation_summary.json',dict(assemblies_prepared=len(prepared),profile_search_executed=False,stage_gate=gate))
            print((out/'preparation_summary.json').read_text());return
        if not gate['allowed']:print(json.dumps(gate,indent=2));return
        run_stage1(cfg,out,args.stage);return
    write_json(out/'stage_gate.json',stage_gate(1,cfg));run_stage1(cfg,out)


if __name__=='__main__':main()
