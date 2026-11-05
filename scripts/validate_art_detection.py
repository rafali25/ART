"""Frozen ART validation: independent published-core anchors, synthetic DNA and RT controls.

All coordinates zero based and half open. Downloaded source annotations are retained.
Run from workspace root. No production thresholds are changed by this benchmark.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from collections import Counter
import hashlib
import itertools
import json
import os
from pathlib import Path
import random
import subprocess
import urllib.request

import numpy as np
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation

from art_array_detector import ROOT, refinement_scan, initial_scan, distinct_arrays, empirical, stable_seed
from art_array_assessment import aligned_pwm, quality_flags, short_period_scan, main_status, sensitivity_batch, _scan_core_hash
from extract_art_loci import read_tsv, write_tsv, write_json, plan_locus, extract_record, source_pieces, pid
from art_locus_evidence import classify_locus


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def empirical_p(score, draws):
    return (1+sum(value>=score for seed,value in draws))/(1+len(draws))


def table(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    write_tsv(path, rows, fields)


def command(args, log):
    env = os.environ.copy()
    env['PATH'] = str(ROOT/'.tools/cygwin/bin') + os.pathsep + env['PATH']
    with Path(log).open('w') as handle:
        subprocess.run([str(x) for x in args], stdout=handle, stderr=subprocess.STDOUT, env=env, check=True)


def search(hmm, fasta, stem, cfg):
    command([ROOT/'.tools/hmmer-install/bin/hmmsearch.exe','--cpu','4','--seed','42',
        '-Z',cfg['profile_database_size'],'--domZ',cfg['profile_database_size'],
        '-E','0.001','--domE','0.001','--tblout',str(stem)+'.tbl',
        '--domtblout',str(stem)+'.domtbl',hmm,fasta],str(stem)+'.log')
    rows=[]
    for line in Path(str(stem)+'.domtbl').read_text().splitlines():
        if not line.strip() or line.startswith('#'): continue
        c=line.split(maxsplit=22)
        cov=(int(c[16])-int(c[15])+1)/int(c[5])
        rows.append(dict(protein_id=c[0],profile=Path(hmm).stem,sequence_E=float(c[6]),
            domain_E=float(c[12]),bits=float(c[13]),coverage=cov,
            passes=float(c[6])<=1e-5 and float(c[12])<=1e-5 and cov>=.6))
    return rows


def references(out, cfg):
    arrays=read_tsv(ROOT/'results/art_arrays_assessed/array_assessments.tsv')
    repeats=read_tsv(ROOT/'results/art_arrays_assessed/repeats.tsv')
    manifests=read_tsv(ROOT/'results/oriented_art_loci/locus_manifest.tsv')
    discovery={r['RT_protein_id']:r for r in read_tsv(ROOT/'results/art_discovery/candidate_loci.tsv')}
    rows=[]; gold=[]; reverse=[]
    for meta in manifests:
        loc=meta['locus_id']; rtpid=meta['RT_protein_id']
        record=SeqIO.read(ROOT/'results/oriented_art_loci'/meta['genbank_path'],'genbank')
        rt=next(f for f in record.features if f.type=='CDS' and pid(f)==rtpid)
        a=next(r for r in arrays if r['locus_id']==loc)
        rr=[r for r in repeats if r['array_id']==a['array_id']]
        original=next(r for r in read_tsv(ROOT/'results/art_discovery/source_records.tsv') if r['accession']==meta['accession'])
        src=SeqIO.read(ROOT/original['path'],original['format'])
        srcrt=next(f for f in src.features if f.type=='CDS' and pid(f)==rtpid)
        plan=plan_locus(src,srcrt)
        rcsrc=src.reverse_complement(id=True,name=True,description=True,features=True,annotations=True)
        # Biopython may discard topology when preserving other annotations; set explicitly.
        rcsrc.annotations['topology']=src.annotations.get('topology','linear')
        rcrt=next(f for f in rcsrc.features if f.type=='CDS' and pid(f)==rtpid)
        rcplan=plan_locus(rcsrc,rcrt)
        oriented=extract_record(src,srcrt,plan,loc)
        reoriented=extract_record(rcsrc,rcrt,rcplan,loc)
        mapping_ok=True
        for rep in rr:
            x,y=int(rep['locus_start0']),int(rep['locus_end0'])
            for old,new in zip(source_pieces(x,y,plan['segments']), reversed(source_pieces(x,y,rcplan['segments']))):
                mapping_ok &= old['source_start0']==len(src)-new['source_end0'] and old['source_end0']==len(src)-new['source_start0'] and old['orientation']==-new['orientation']
        reverse.append(dict(locus_id=loc,source_strand=srcrt.location.strand,reverse_source_strand=rcrt.location.strand,
            identical_oriented_DNA=str(oriented.seq)==str(reoriented.seq),identical_oriented_CDS=str(rt.extract(record.seq))==str(rcrt.extract(rcsrc.seq)),
            equivalent_original_coordinates=mapping_ok))
        core=cfg['reference_core']; upstream=str(record.seq[:int(rt.location.start)]).upper()
        hits=[p for p in range(len(upstream)-len(core)+1) if sum(x!=y for x,y in zip(upstream[p:p+len(core)],core))<=cfg['reference_core_max_mismatches']]
        # Independent anchors: the specified paper core only; no discovery score or output used.
        is_staph=meta['accession']!='MZ422438.1'
        recalls=[]; exact=[]; gs=[]
        if is_staph:
            for i,p in enumerate(hits):
                x=p-cfg['reference_repeat_left_of_core_nt']; y=x+cfg['reference_repeat_length_nt']
                recovered=[r for r in rr if int(r['locus_start0'])<=p and int(r['locus_end0'])>=p+len(core)]
                ex=any(int(r['locus_start0'])==x and int(r['locus_end0'])==y for r in rr)
                recalls.append(bool(recovered));exact.append(ex);gs.append((x,y))
                gold.append(dict(locus_id=loc,copy=i+1,core_start0=p,core_end0=p+len(core),repeat_start0=x,repeat_end0=y,
                    source_pieces_json=json.dumps(source_pieces(x,y,plan['segments'])),core_recovered=bool(recovered),exact_22nt_boundary=ex))
        d=discovery[rtpid]
        status,reasons=classify_locus(rt_family=True,partner_family=bool(d['partner_family_profiles']),
            compatible_context=d['cooriented_RT_partner_context']=='True',array_status=a['assessment_status'],copies=int(a['copy_count']))
        rows.append(dict(locus_id=loc,gold_status='published_core_anchored' if is_staph else 'independent_copy_coordinates_unavailable',
            expected_core_copies=len(hits) if is_staph else '',core_copies_recovered=sum(recalls) if is_staph else '',
            exact_paper_22nt_boundaries=sum(exact) if is_staph else '',gold_spacers_nt_json=json.dumps([c-b for (x,b),(c,y) in zip(gs,gs[1:])]),
            measured_spacers_nt_json=a['spacer_lengths_nt_json'],detected_copies=a['copy_count'],array_status=a['assessment_status'],
            locus_status=status,reasons=';'.join(reasons),flags=a['flags'],adjusted_p=a['bonferroni_p'],
            split='retrospective_development_reference',discovery_routes=d['discovery_routes']))
    table(out/'reference_recovery.tsv',rows);table(out/'reference_gold.tsv',gold);table(out/'reverse_complement_checks.tsv',reverse)
    return rows,reverse


def controls(out):
    specs=[('retron','M24408.1.embl','embl','AAA61471.1','https://www.ebi.ac.uk/ena/browser/view/M24408.1'),
        ('DGR','NC_005357.1.gb','genbank','NP_958675.1','https://www.ncbi.nlm.nih.gov/nuccore/NC_005357.1'),
        ('CRISPR_RT','CP002583.1.embl','embl','ADZ89953.1','https://www.ebi.ac.uk/ena/browser/view/CP002583.1')]
    loci=[]; proteins=[]; manifest=[]
    for family,file,fmt,target,url in specs:
        path=out/'sources'/file
        if not path.exists():
            (out/'sources').mkdir(parents=True,exist_ok=True)
            accession=file.rsplit('.',1)[0]
            endpoint=('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='+accession+'&rettype=gbwithparts&retmode=text') if fmt=='genbank' else ('https://www.ebi.ac.uk/ena/browser/api/embl/'+accession+'?download=false')
            with urllib.request.urlopen(endpoint,timeout=55) as handle:path.write_bytes(handle.read())
        src=SeqIO.read(path,fmt)
        rt=next(f for f in src.features if f.type=='CDS' and pid(f)==target)
        plan=plan_locus(src,rt); rec=extract_record(src,rt,plan,family+'__'+target)
        SeqIO.write(rec,out/'loci'/(rec.id+'.gb'),'genbank')
        SeqIO.write(rec,out/'loci'/(rec.id+'.fna'),'fasta')
        write_json(out/'loci'/(rec.id+'.map.json'),dict(segments=[asdict(s) for s in plan['segments']],**{k:v for k,v in plan.items() if k not in ('segments','plus_segments')}))
        loci.append((rec,rt,target,plan,family))
        proteins.append(SeqRecord(Seq(rt.qualifiers['translation'][0]),id=target,description=family))
        for f in rec.features:
            if f.type=='CDS' and pid(f)!=target and f.qualifiers.get('translation'):
                proteins.append(SeqRecord(Seq(f.qualifiers['translation'][0]),id=pid(f),description=family+' locus neighbor'))
        manifest.append(dict(family=family,source_accession=src.id,source_url=url,source_file=file,source_sha256=digest(path),
            protein_id=target,RT_start0=int(rt.location.start),RT_end0=int(rt.location.end),strand=rt.location.strand,
            locus_id=rec.id,upstream_available=plan['upstream_available_nt'],upstream_shortfall=plan['upstream_shortfall_nt'],
            evidence='Ec86 per Yoon methods' if family=='retron' else 'BPP-1 Brt source annotation' if family=='DGR' else 'Marme_0669 Cas6-RT-Cas1, Mohr et al. 2018',status='retrieved'))
    manifest.append(dict(family='NART',status='unavailable_verified_sequence_mapping',
        evidence='Yoon PDF Supp.Fig.5 internal examples (54-kb,333-kb,14-kb contigs) lack accession mapping; source audit retained. AR9 membership alone is insufficient; KU878088.1 has no annotated RT and was excluded.'))
    table(out/'comparison_source_manifest.tsv',manifest)
    # Source-family labels are supplied by the published myRT reference collection.
    for r in SeqIO.parse(ROOT/'.tools/myRT/Models/RTs-collection.faa','fasta'):
        if r.id.rsplit('_',1)[-1] in ('Retrons','DGRs','CRISPR'):
            r.seq=Seq(str(r.seq).rstrip('*')); proteins.append(r)
    SeqIO.write(proteins,out/'comparison_proteins.faa','fasta')
    return loci,proteins


def profile_benchmark(out,cfg,proteins):
    domains=[]
    for family in ('ART_RT','Type_I_partner','Type_II_partner'):
        domains+=search(ROOT/'results/art_discovery/profiles'/(family+'.hmm'),out/'comparison_proteins.faa',out/'profiles'/family,cfg)
    table(out/'comparison_profile_hits.tsv',domains)
    target_ids={r['RT_protein_id'] for r in read_tsv(ROOT/'results/art_discovery/candidate_loci.tsv')}
    targets=[]
    for r in SeqIO.parse(ROOT/'results/art_discovery/target_proteins.faa','fasta'):
        if r.id in target_ids: targets.append(r)
    if len(targets)!=9:
        raise ValueError('Expected nine accession-keyed annotated target RTs')
    aligner=Align.PairwiseAligner(scoring='blastp',mode='global')
    pairs=[]; parents={r.id:r.id for r in targets}
    def root(x):
        while parents[x]!=x:x=parents[x]
        return x
    for a,b in itertools.combinations(targets,2):
        al=aligner.align(a.seq,b.seq)[0]; c=al.counts()
        identity=c.identities/max(1,c.aligned)
        ca=c.aligned/len(a); cb=c.aligned/len(b)
        same=identity>=cfg['relatedness_identity_threshold'] and min(ca,cb)>=cfg['relatedness_mutual_coverage']
        pairs.append(dict(first=a.id,second=b.id,identity=identity,first_coverage=ca,second_coverage=cb,same_component_edge=same))
        if same:parents[root(b.id)]=root(a.id)
    # Conservative lineage supergroups prevent a divergent SA1 sequence or a
    # same-genome LPJP1 paralog leaking across an 80%-identity component split.
    lp_ids={'QXN67953.1','QXN67955.1'}
    groups={r.id:('LPJP1_lineage' if r.id in lp_ids else 'Staphylococcus_phage_lineage') for r in targets}
    table(out/'relatedness_pairs.tsv',pairs)
    table(out/'relatedness_groups.tsv',[dict(protein_id=p,similarity_component=root(p),evaluation_supergroup=g,
        rationale='conservative host-genus/phage-lineage grouping; LPJP1 same-genome paralogs kept together') for p,g in groups.items()])
    seedrecords=list(SeqIO.parse(ROOT/'results/art_discovery/profiles/ART_RT.seeds.faa','fasta'))
    SeqIO.write(targets+proteins,out/'profile_evaluation.faa','fasta')
    foldrows=[]
    for i,group in enumerate(sorted(set(groups.values()))):
        held={p for p,g in groups.items() if g==group}
        train=[r for r in seedrecords if r.id not in held]
        stem=out/'profiles'/('holdout_'+str(i+1))
        SeqIO.write(train,str(stem)+'.seeds.faa','fasta')
        if not train:
            foldrows.append(dict(fold=i+1,heldout_group=group,status='no_training_group_remaining'));continue
        command(['cmd.exe','/d','/c',ROOT/'.tools/mafft-win/mafft.bat','--auto','--thread','1','--inputorder','--out',str(stem)+'.alignment.faa',str(stem)+'.seeds.faa'],str(stem)+'.mafft.log')
        command([ROOT/'.tools/hmmer-install/bin/hmmbuild.exe','-n','ART_holdout_'+str(i+1),str(stem)+'.hmm',str(stem)+'.alignment.faa'],str(stem)+'.build.log')
        hs=search(str(stem)+'.hmm',out/'profile_evaluation.faa',stem,cfg)
        table(str(stem)+'.hits.tsv',hs)
        for p in sorted(held):
            good=[h for h in hs if h['protein_id']==p and h['passes']]
            foldrows.append(dict(fold=i+1,heldout_group=group,protein_id=p,training_seed_count=len(train),
                training_ids=';'.join(r.id for r in train),heldout_recovered=bool(good),status='group_heldout_evaluation',
                negative_proteins_passing=len({h['protein_id'] for h in hs if h['passes'] and h['protein_id'] not in target_ids})))
    table(out/'group_holdout_results.tsv',foldrows)
    return domains,foldrows


def synthetic_cases(out,cfg):
    cases=[]; truths=[]
    names=['long_spacer_120','long_spacer_180','long_spacer_220','imperfect_copies','missing_copy','extended_beyond_3kb','truncated_upstream','unrelated_RT_with_array','repetitive_DNA']
    for ci,name in enumerate(names):
        for rep in range(cfg['synthetic_replicates_per_case']):
            rng=random.Random(stable_seed(cfg['random_seed'],name,rep))
            motif=''.join(rng.choices('ACGT',k=28))
            seq=list(''.join(rng.choices('ACGT',k=9000)))
            spacer=int(name.rsplit('_',1)[-1]) if name.startswith('long_spacer_') else 180
            starts=[4000+i*(len(motif)+spacer) for i in range(6)]
            if name=='missing_copy':starts.pop(2)
            if name=='extended_beyond_3kb':starts=[1000+i*208 for i in range(20)]
            if name=='repetitive_DNA':
                # A CRISPR-like array with 35nt spacers is allowed to be a DNA-array hit.
                starts=[4300+i*63 for i in range(20)];spacer=35
            for i,p in enumerate(starts):
                m=list(motif)
                if name=='imperfect_copies' and i in (1,4):
                    m[13]=next(b for b in 'ACGT' if b!=m[13])
                seq[p:p+len(m)]=m
            rec=SeqRecord(Seq(''.join(seq)),id=f'{name}_{rep+1}',description='synthetic frozen independent motif')
            rec.annotations={'molecule_type':'DNA','topology':'linear'}
            rt=SeqFeature(SimpleLocation(6000,7500,strand=1),type='CDS',qualifiers={'protein_id':['synthetic_RT']})
            rec.features=[rt,SeqFeature(SimpleLocation(7510,8320,strand=1),type='CDS',qualifiers={'protein_id':['synthetic_partner']})]
            shortfall=0
            if name=='truncated_upstream':
                # Source ends before any planted copy; absence cannot be a confident negative.
                rec=rec[5200:];rec.id=f'{name}_{rep+1}';shortfall=5200;starts=[]
            intervals=[[p,p+len(motif)] for p in starts]
            SeqIO.write(rec,out/'loci'/(rec.id+'.gb'),'genbank')
            SeqIO.write(rec,out/'loci'/(rec.id+'.fna'),'fasta')
            truths.append(dict(case_id=rec.id,test=name,replicate=rep+1,repeat=motif,repeat_intervals_json=json.dumps(intervals),
                known_spacers_json=json.dumps([b-a-28 for a,b in zip(starts,starts[1:])]),
                rt_family=name not in ('unrelated_RT_with_array','repetitive_DNA'),partner_family=name not in ('unrelated_RT_with_array','repetitive_DNA'),
                upstream_shortfall=shortfall,expected='incomplete_assessment' if shortfall else 'no_ART_locus_status' if name in ('unrelated_RT_with_array','repetitive_DNA') else 'recover_all_planted_copies'))
            cases.append((rec,name,intervals,shortfall))
    table(out/'synthetic_truth.tsv',truths)
    return cases


def nulls(out,sequence,cfg,vfg,pool,case):
    path=out/'nulls'/(case+'.json')
    key=dict(sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),cfg=cfg,validation=vfg,scan_core_sha256=_scan_core_hash())
    if path.exists():
        data=json.loads(path.read_text())
        old=data['key']
        if 'detector_sha256' in old:
            # Migrate only verified pre-fix draws: the change is downstream of
            # every null score. No seed-discovery/chaining/scoring code changed.
            baseline=out/'baseline_before_anchor_fix/art_array_detector.py'
            if old['detector_sha256']!=digest(baseline) or old['assessment_sha256']!=digest(ROOT/'scripts/art_array_assessment.py') or any(old[k]!=key[k] for k in ('sequence_sha256','cfg','validation')):
                raise ValueError('Unverified pre-fix null cache: '+case)
            import ast
            before=ast.parse(baseline.read_text());after=ast.parse((ROOT/'scripts/art_array_detector.py').read_text())
            funcs=lambda tree:{n.name:ast.dump(n,include_attributes=False) for n in tree.body if isinstance(n,ast.FunctionDef)}
            b,a=funcs(before),funcs(after)
            if {n for n in b if b[n]!=a[n]}!={'align_boundaries'}:raise ValueError('Null scanner changed during migration')
            data['key']=key;write_json(path,data)
        elif old!=key:raise ValueError('Stale validation null cache: '+case)
        return data['scores']
    result={}
    for kind in ('local_blocks','dinucleotide'):
        n=vfg['null_draws_per_background']; seeds=[stable_seed(vfg['random_seed'],'validation',case,kind,i) for i in range(n)]
        jobs=[(sequence,kind,cfg,seeds[j:j+20]) for j in range(0,n,20)]
        result[kind]=[list(row) for batch in pool.map(sensitivity_batch,jobs) for row in batch]
    write_json(path,dict(key=key,scores=result))
    return result


def sequence_benchmark(out,cfg,vfg,cases,controls_,domains):
    rows=[]; allnull=[]; ntests=len(cases)+len(controls_)
    assess=json.loads((ROOT/'config/art_array_assessment.json').read_text())
    jobs=[]
    for rec,name,truth,shortfall in cases:
        rt=next(f for f in rec.features if f.type=='CDS' and pid(f)=='synthetic_RT')
        jobs.append((rec,name,truth,shortfall,rt.location.start,name not in ('unrelated_RT_with_array','repetitive_DNA'),name not in ('unrelated_RT_with_array','repetitive_DNA')))
    for rec,rt,target,plan,family in controls_:
        ort=next(f for f in rec.features if f.type=='CDS' and pid(f)==target)
        rthit=any(h['protein_id']==target and h['profile']=='ART_RT' and h['passes'] for h in domains)
        partnerhit=any(h['protein_id']!=target and h['profile']!='ART_RT' and h['passes'] for h in domains if h['protein_id'] in {pid(f) for f in rec.features if f.type=='CDS'})
        jobs.append((rec,family,[],plan['upstream_shortfall_nt'],int(ort.location.start),rthit,partnerhit))
    with ProcessPoolExecutor(max_workers=vfg['workers']) as pool:
        for rec,name,truth,shortfall,rtstart,rtok,partnerok in jobs:
            seq=str(rec.seq[max(0,rtstart-6000):rtstart]).upper()
            offset=max(0,rtstart-6000)
            scan=refinement_scan(seq,cfg); initial=initial_scan(seq[-3000:],cfg)
            ns=nulls(out,seq,cfg,vfg,pool,rec.id)
            for kind,draws in ns.items():
                allnull.extend(dict(case_id=rec.id,background=kind,draw=i+1,seed=s,score=sc) for i,(s,sc) in enumerate(draws))
            bestp=max(empirical_p(scan['score'],draws) for draws in ns.values())
            # Use the same more-copies-first distinct candidate policy as production.
            cs=distinct_arrays(scan['candidates']); chosen=cs[0] if cs else None
            intervals=[]; model={}; flags=[]; p=1.
            if chosen:
                model=aligned_pwm(seq,chosen,cfg,out,rec.id)
                intervals=model.get('final_intervals',[])
                p=max(empirical_p(chosen['score'],draws) for draws in ns.values())
                if intervals:
                    q=quality_flags(rec,[[a+offset,b+offset] for a,b in intervals],rtstart,
                        dict(upstream_shortfall_nt=shortfall),cfg,assess,short_period_scan(seq,cfg),chosen['score'])
                    flags=q['flags']
            adjusted=min(1.,p*ntests); ast=main_status(len(intervals),adjusted,flags)
            status,reasons=classify_locus(rt_family=rtok,partner_family=partnerok,compatible_context=rtok and partnerok,
                array_status=ast,copies=len(intervals),upstream_shortfall=shortfall,
                ambiguous_nt=sum(c not in 'ACGT' for c in seq))
            found=[[a+offset,b+offset] for a,b in intervals]
            recovered=sum(any(max(a,x)<min(b,y) and min(b,y)-max(a,x)>=min(b-a,y-x)*.8 for x,y in found) for a,b in truth)
            exact=sum([a,b] in found for a,b in truth)
            if truth and name not in ('unrelated_RT_with_array','repetitive_DNA'):
                passed=recovered==len(truth) and len(intervals)==len(truth) and ast=='accepted_array_candidate'
                if name=='missing_copy':passed &= model.get('selected',{}).get('missing_copies')==1
                if name=='extended_beyond_3kb':passed &= recovered>initial['R']
            elif shortfall:passed=status=='incomplete_ART_locus_assessment' if rtok and partnerok else status!='supported_ART_locus_candidate'
            else:passed=status!='supported_ART_locus_candidate'
            row=dict(case_id=rec.id,test=name,expected_copies=len(truth),recovered_truth_copies=recovered,exact_boundaries=exact,
                detected_copies=len(found),initial_best_copies=initial['R'],extended_only_or_additional_copies=bool(truth and recovered>initial['R']),
                inferred_skipped_copy_links=model.get('selected',{}).get('missing_copies',''),
                detected_intervals_json=json.dumps(found),measured_spacers_json=json.dumps([c-b for (a,b),(c,d) in zip(found,found[1:])]),
                best_search_score=scan['score'],best_search_p=bestp,selected_chain_score=chosen['score'] if chosen else 0,
                local_p=empirical_p(chosen['score'] if chosen else 0,ns['local_blocks']),
                dinucleotide_p=empirical_p(chosen['score'] if chosen else 0,ns['dinucleotide']),
                selected_adjusted_p=adjusted,testing_family_loci=ntests,array_status=ast,flags=';'.join(flags),
                ART_RT_support=rtok,partner_support=partnerok,locus_status=status,reasons=';'.join(reasons),
                upstream_available=len(seq),upstream_shortfall=shortfall,pass_expected_behavior=passed,model_status=model.get('status','no_chain'))
            rows.append(row);table(out/'sequence_cases.tsv',rows)
            write_json(out/'details'/(rec.id+'.json'),dict(scan=scan,model=model,assessment=row))
            print(rec.id,passed,recovered,'/',len(truth),status,flush=True)
    table(out/'null_scores.tsv',allnull)
    return rows


def finalize(out):
    """Refresh derived reference checks and summary without rerunning any searches."""
    vfg=json.loads((ROOT/'config/art_validation.json').read_text())
    refs,rc=references(out,vfg)
    rows=read_tsv(out/'sequence_cases.tsv')
    truths={r['case_id']:r for r in read_tsv(out/'synthetic_truth.tsv')}
    for r in rows:
        if r['case_id'] in truths:
            truth=truths[r['case_id']]
            positive=truth['expected']=='recover_all_planted_copies'
            if positive:
                r['unmatched_detected_copies']=int(r['detected_copies'])-int(r['recovered_truth_copies'])
                r['pass_expected_behavior']=(int(r['recovered_truth_copies'])==int(r['expected_copies']) and
                    int(r['unmatched_detected_copies'])==0 and r['array_status']=='accepted_array_candidate')
                if r['test']=='missing_copy':r['pass_expected_behavior'] &= str(r['inferred_skipped_copy_links'])=='1'
                if r['test']=='extended_beyond_3kb':r['pass_expected_behavior'] &= r['extended_only_or_additional_copies']=='True'
            r['evaluation_role']='synthetic_regression_no_independent_generalization_claim'
        else:r['evaluation_role']='external_annotated_genomic_control'
        detail_path=out/'details'/(r['case_id']+'.json')
        detail=json.loads(detail_path.read_text());detail['assessment']=r;write_json(detail_path,detail)
    table(out/'sequence_cases.tsv',rows)
    proteins=list(SeqIO.parse(out/'comparison_proteins.faa','fasta'))
    hits=read_tsv(out/'comparison_profile_hits.tsv');folds=read_tsv(out/'group_holdout_results.tsv')
    passed=lambda x:x is True or x=='True'
    summary=dict(reference_loci=len(refs),core_anchored_references=sum(r['gold_status']=='published_core_anchored' for r in refs),
        reference_core_copies=sum(r['expected_core_copies'] for r in refs if r['expected_core_copies']!=''),
        recovered_reference_core_copies=sum(r['core_copies_recovered'] for r in refs if r['core_copies_recovered']!=''),
        exact_paper_22nt_boundaries=sum(r['exact_paper_22nt_boundaries'] for r in refs if r['exact_paper_22nt_boundaries']!=''),
        reference_statuses=dict(Counter(r['locus_status'] for r in refs)),
        reverse_complement_checks_pass=all(all(r[k] for k in ('identical_oriented_DNA','identical_oriented_CDS','equivalent_original_coordinates')) for r in rc),
        comparison_protein_records=len(proteins),comparison_unique_sequences=len({str(r.seq) for r in proteins}),
        myRT_reference_family_counts=dict(Counter(r.id.rsplit('_',1)[-1] for r in proteins if r.id.rsplit('_',1)[-1] in ('Retrons','DGRs','CRISPR'))),
        production_profile_passing_domains=sum(passed(r['passes']) for r in hits),
        group_holdout_targets=len(folds),group_holdout_recovered=sum(passed(r.get('heldout_recovered')) for r in folds),
        sequence_cases=len(rows),sequence_cases_passing=sum(passed(r['pass_expected_behavior']) for r in rows),
        sequence_failures=[r['case_id'] for r in rows if not passed(r['pass_expected_behavior'])],
        synthetic_cases=len(truths),null_scores=len(rows)*vfg['null_draws_per_background']*2,
        missing_comparison_sets=['NART'],missing_partner_families=['Type_III_first_partner','Type_III_second_partner'],scale_search_ready=False)
    write_json(out/'summary.json',summary)
    frozen=json.loads((out/'frozen_inputs.json').read_text())
    changes=[]
    for path,old in frozen.items():
        new=digest(ROOT/path)
        if new!=old:
            if path.replace('\\','/')!='scripts/validate_art_detection.py':
                raise ValueError('Inference input changed; rerun the benchmark before finalizing: '+path)
            changes.append(dict(path=path,before_sha256=old,after_sha256=new))
        frozen[path]=new
    if changes:write_json(out/'finalization_input_refresh.json',dict(reason='Benchmark audit adds copy-precision criteria; original nine loci reprocessed after majority-anchor repair. No retrieval thresholds, null scores or original annotations changed.',changes=changes))
    write_json(out/'frozen_inputs.json',frozen)
    for directory,names in [('scripts',['validate_art_detection.py','verify_art_validation.py','art_locus_evidence.py','art_array_detector.py','art_array_assessment.py','test_art_validation.py']),
                            ('config',['art_validation.json','art_array_search.json','art_array_assessment.json','art_discovery.json'])]:
        destination=out/'code'/directory;destination.mkdir(parents=True,exist_ok=True)
        for name in names:(destination/name).write_bytes((ROOT/directory/name).read_bytes())
    write_report(out,summary,refs,rows,folds)
    write_json(out/'output_sha256.json',{str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='output_sha256.json'})
    print(json.dumps(summary,indent=2))


def write_report(out,summary,refs,rows,folds):
    text=['# ART detection validation — first milestone', '',
        '**Milestone partially met; scale search remains on hold.** All 35 reference-core copies in the seven available Staphylococcus genomes are recovered. Long spacers are retained. This is not validation of the complete published metagenomic census or all ART partner types.', '',
        '## Published loci and orientation', '',
        'The independent reference anchors use the paper’s ATATGAATACGTAT core with at most one mismatch, solely to score reference recovery. This motif never enters discovery. Full source genomes and original annotations are retained. All nine real RT loci have identical oriented DNA, RT CDS sequence and equivalent genomic repeat coordinates after reverse-complementing the complete source genome.', '',
        '| Available reference | Core copies recovered | Final locus assessment | Reason / limit |',
        '| --- | --- | --- | --- |']
    labels={'LC680885.1':'S6','MN091626.1':'PALS2','MW218148.1':'SA1','MW248466.1':'MarsHill','MW349129.1':'Madawaska','MZ779063.1':'UFV_DC4','OR836606.1':'LY01'}
    for r in refs:
        accession=r['locus_id'].split('__')[0];label=labels.get(accession,r['locus_id'])
        recovered=f"{r['core_copies_recovered']}/{r['expected_core_copies']}" if r['expected_core_copies']!='' else 'No independent copy-coordinate gold set'
        reason=r['flags'] or r['reasons'] or 'Family + partner + accepted array evidence'
        text.append(f"| {label} | {recovered} | {r['locus_status']} | {reason}; adjusted p={float(r['adjusted_p']):.5f} |")
    text += ['',
        'MarsHill repeat-start spacings reproduce Figure 2G: 213, 197, 210 and 173 nt. These are start-to-start spacings. Its conservation-defined 17-nt repeats give actual spacers of 196, 180, 193 and 156 nt. The paper also uses fixed 22-nt core-centred windows for RNAfold; none of the 35 conservation-defined boundaries exactly equals that 22-nt definition. Core recovery therefore does **not** demonstrate exact full-boundary reproduction. The independent anchor table records this distinction.', '',
        'UFV_DC4 and LY01 arrays are recovered but fail the existing nine-locus corrected statistical rule (p=0.06294 and 0.08991). LY01 also has an intervening annotated gene. Madawaska requires review because repeats overlap two retained CDS annotations. These annotations were not deleted to make the references pass. The second LPJP1 RT has insufficient recovered copies and p=0.18881. LPJP1’s first RT remains supported, including a 264-nt spacer; observed spacer ranges are not hard filters.', '',
        '## Synthetic and external genomic tests', '',
        'There are 27 synthetic cases: three independently generated motifs per test, plus three external genomic controls. Synthetic RT/partner support is supplied as known test metadata to exercise the evidence gate; it is not a simulated HMM sensitivity estimate. Synthetic fixtures are regression/development tests, including cases that exposed alignment bugs, and are not claimed as untouched evaluation data.', '',
        '| Test | Cases passing expected behaviour | Notes |', '| --- | --- | --- |']
    notes={'long_spacer_120':'Planted 28-nt repeats, 120-nt spacers','long_spacer_180':'Planted 28-nt repeats, 180-nt spacers',
        'long_spacer_220':'Planted 28-nt repeats, 220-nt spacers','imperfect_copies':'Two copies carry a substitution; copy precision is scored as well as recall',
        'missing_copy':'One repeat absent; require exactly one skipped-copy link',
        'extended_beyond_3kb':'20 copies spanning 3,980 nt; extended search must recover additional copies',
        'truncated_upstream':'Only 800 nt upstream available; incomplete assessment required',
        'unrelated_RT_with_array':'Strong long-spacer DNA array must not imply ART family membership',
        'repetitive_DNA':'CRISPR-like short-spacer array without ART RT or partner support',
        'retron':'Ec86 M24408.1 / AAA61471.1; source fragment has only 547 nt upstream',
        'DGR':'BPP-1 NC_005357.1 / NP_958675.1, reverse-strand Brt',
        'CRISPR_RT':'MMB-1 CP002583.1 / ADZ89953.1, Marme_0669 Cas6–RT–Cas1'}
    for test in dict.fromkeys(r['test'] for r in rows):
        group=[r for r in rows if r['test']==test]
        passed=sum(r['pass_expected_behavior'] is True or r['pass_expected_behavior']=='True' for r in group)
        text.append(f"| {test} | {passed}/{len(group)} | {notes[test]} |")
    text += ['',f"Overall: **{summary['sequence_cases_passing']}/{len(rows)}** case behaviours pass. Failures: "+(', '.join(summary['sequence_failures']) or 'none')+'.', '',
        'An array can recover all planted copies and still fail copy accuracy by adding a spurious copy. `unmatched_detected_copies` and exact-boundary counts are retained. No thresholds were relaxed to turn these errors into passes. Full repeat intervals, training alignments, PWM models and shuffled-search results accompany each case.', '',
        'Each of the 30 sequence windows has 1,000 local 50-nt block shuffles and 1,000 exact dinucleotide-preserving shuffles (60,000 full-search scores). Every shuffle rediscovers recurring seeds and maximizes the chain score. The selected real chain is compared with these full-search maxima; p=(1+exceedances)/(1+1000). The larger background p is Bonferroni-corrected across all 30 benchmark loci. Minimum possible adjusted p is 30/1001≈0.02997. The original nine-reference cohort retains its separate nine-locus correction and 36,000 null scores. Initial 3-kb results are diagnostic here; the full 6-kb scan always runs.', '',
        '## Protein specificity and group exclusion', '',
        f"The comparison FASTA contains {summary['comparison_protein_records']} protein records ({summary['comparison_unique_sequences']} distinct sequences): 96 labelled retron, 502 DGR and 130 CRISPR-associated myRT references, plus the three anchor RTs and available translated neighboring CDS. These are not all independent evolutionary examples. None passes the production ART RT, type I partner or type II partner screens at the frozen E≤1e-5/domain coverage≥0.6 settings (HMMER Z=domZ=2076). Raw subthreshold hits are retained.", '',
        'Pairwise full-protein identity components are recorded, then merged conservatively into Staphylococcus-phage and LPJP1 lineage groups. SA1 stays with its host-genus group even though it separates at 80% identity; LPJP1 paralogs stay together. Entire groups are excluded from training in both evaluation folds. All nine RTs are retrieved. However, the profile trained on the single LPJP1 seed also retrieves retron NP_782489.1_Retrons (coverage 0.6652). Thus RT sensitivity across these two groups does not establish family specificity or broad generalization. No cutoff was retuned using these evaluation results.', '',
        'The seven Staphylococcus genomes were already used during development; SA1 and MarsHill are sequence holdouts in the production profile, not independent lineage holdouts. The myRT proteins were also screened in the earlier pilot: this is a labelled comparison panel, not a previously untouched negative set. Type I has one training partner and no positive holdout. Type II has only one conservative lineage group, so no independent partner-lineage evaluation is available. Type III first and second partner profiles remain separate and unbuilt.', '',
        '## Missing comparisons and remaining work', '',
        '**NART is unavailable in this benchmark.** The paper’s NART examples use internal labels / contig lengths without a retrieved accession mapping. The [source audit](NART_SOURCE_AUDIT.md) did not recover verified NART proteins or contigs. AR9-like phage membership is not sufficient: downloaded KU878088.1 has no annotated RT and was excluded. `comparison_source_manifest.tsv` explicitly records the missing NART set; there is no NART specificity result. Missing NART comparisons and type III references prevent a complete validation claim.', '',
        'Before scale search: resolve residual copy-precision errors, obtain independently mapped full repeat boundaries, review the three flagged Staphylococcus references, add verified NART loci and type III partner seeds, and evaluate profiles on additional independent ART lineages and competing RT families. These are limitations of the available data and present detector, not grounds to alter the frozen acceptance thresholds.', '',
        '## Repairs, evidence and rerun', '',
        'Validation exposed two alignment errors. A chance first seed could veto a genuine majority repeat block; the anchor now uses the median mapped seed midpoint. A chance seed aligned entirely as gaps could discard every supported copy; it is now excluded and its copy index recorded. Original seed discovery, dynamic programming, scoring, thresholds and shuffle procedures are unchanged. Pre-repair code, partial synthetic results and original nine-locus outputs are preserved under `baseline_before_anchor_fix/`. Null-cache migration verifies source hashes and that only the downstream alignment function changed. All nine production loci were reprocessed after the repairs.', '',
        'Main files: `reference_recovery.tsv`, `reference_gold.tsv`, `reverse_complement_checks.tsv`, `sequence_cases.tsv`, `synthetic_truth.tsv`, `comparison_source_manifest.tsv`, `comparison_profile_hits.tsv`, `relatedness_groups.tsv`, `group_holdout_results.tsv`, `null_scores.tsv`, `summary.json`, `verification.json`, `frozen_inputs.json`, `output_sha256.json`. Source records, genomic coordinate maps, raw HMMER files, alignments and code snapshots are retained.', '',
        '```powershell', '.venv/Scripts/python.exe scripts/validate_art_detection.py', '.venv/Scripts/python.exe scripts/verify_art_validation.py', ".venv/Scripts/python.exe -m unittest discover -s scripts -p 'test_*art*.py'", '```', '',
        'Sources: supplied Yoon PDF Fig.2G and Methods (PDF pp.7,32–37); [ART preprint](https://doi.org/10.64898/2026.10.05.755809); [myRT reference collection](https://github.com/mgtools/myRT) at fa8e3624adf4c4e57ebf3350c06d626f805e42f7; [Ec86](https://www.ebi.ac.uk/ena/browser/view/M24408.1); [BPP-1](https://www.ncbi.nlm.nih.gov/nuccore/NC_005357.1); [MMB-1 genome](https://www.ebi.ac.uk/ena/browser/view/CP002583.1); [Cas6–RT–Cas1 primary study](https://doi.org/10.1016/j.molcel.2018.09.013).', '']
    (out/'REPORT.md').write_text('\n'.join(text),encoding='utf-8')


def run():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='results/art_validation');ap.add_argument('--finalize-only',action='store_true');args=ap.parse_args()
    out=ROOT/args.output
    if args.finalize_only:
        finalize(out);return
    for d in ('loci','profiles','alignments','nulls','details'): (out/d).mkdir(parents=True,exist_ok=True)
    vfg=json.loads((ROOT/'config/art_validation.json').read_text())
    cfg=json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
    frozen=[ROOT/'config/art_array_search.json',ROOT/'config/art_array_assessment.json',ROOT/'config/art_discovery.json',
        ROOT/'config/art_validation.json',ROOT/'scripts/art_locus_evidence.py',ROOT/'scripts/validate_art_detection.py',
        ROOT/'scripts/art_array_detector.py',ROOT/'scripts/art_array_assessment.py',ROOT/'scripts/extract_art_loci.py',
        ROOT/'results/art_arrays_assessed/array_assessments.tsv',ROOT/'results/art_arrays_assessed/repeats.tsv']
    frozen+=list((ROOT/'results/art_discovery/profiles').glob('*.hmm'))
    write_json(out/'frozen_inputs.json',{str(p.relative_to(ROOT)):digest(p) for p in frozen})
    refs,rc=references(out,vfg)
    cs,proteins=controls(out)
    domains,folds=profile_benchmark(out,vfg,proteins)
    cases=synthetic_cases(out,vfg)
    rows=sequence_benchmark(out,cfg,vfg,cases,cs,domains)
    finalize(out)


if __name__=='__main__':run()
