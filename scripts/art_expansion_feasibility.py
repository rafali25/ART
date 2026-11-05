"""Separate sequence trees and a transparent necessary-condition block assessment."""
from collections import Counter, defaultdict
import hashlib
import itertools
import json
import re
from pathlib import Path
import time

from Bio import Align, AlignIO, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from art_expansion_workflow import ROOT, OUT, OLD, CPUS, PROFILES, sha, table, read_tsv, write_json, run, checkpoint
from art_reconstruction_common import parse_domains
from extract_art_loci import Segment,source_pieces


def source_hash(path,fmt,expected=None):
    if expected:assert sha(ROOT/str(path).replace('\\','/'))==expected,'Source file hash changed'
    record=SeqIO.read(ROOT/str(path).replace('\\','/'),fmt)
    dna=str(record.seq).upper()
    return min(hashlib.sha256(dna.encode()).hexdigest(),hashlib.sha256(str(Seq(dna).reverse_complement()).encode()).hexdigest())


def canonical_RT_interval(r):
    source=SeqIO.read(ROOT/r['source_path'].replace('\\','/'),r['source_format'])
    locus=SeqIO.read(ROOT/r['locus_genbank'].replace('\\','/'),'genbank')
    mapping=json.loads((ROOT/r['coordinate_map'].replace('\\','/')).read_text())
    segments=[Segment(**s) for s in mapping['segments']]
    feature=next(f for f in locus.features if f.type=='CDS' and f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',['']))[0]==r['RT_feature_id'])
    parts=[s for p in feature.location.parts for s in source_pieces(int(p.start),int(p.end),segments)]
    if len(parts)!=1:return None
    a,b,strand=parts[0]['source_start0'],parts[0]['source_end0'],parts[0]['orientation']
    if hashlib.sha256(str(source.seq).upper().encode()).hexdigest()!=r['source_DNA_sha256']:
        a,b,strand=len(source)-b,len(source)-a,-strand
    return a,b,strand


def collect_pairs():
    result=[]
    base=ROOT/'results/art_type1_catalog'
    rs={r.id:str(r.seq) for r in SeqIO.parse(base/'RT_proteins.faa','fasta')}
    ps={r.id:str(r.seq) for r in SeqIO.parse(base/'primary_partners.faa','fasta')}
    for r in read_tsv(base/'typeI_samples.tsv'):
        lid=r['locus_id']
        cohort={v['locus_id']:v for v in read_tsv(ROOT/r['result_folder']/'candidate_loci.tsv')}
        source_meta={v['accession']:v for v in read_tsv(ROOT/r['result_folder']/'source_records.tsv')}
        source=source_meta[r['contig_id']]
        profile=json.loads(cohort[lid].get('art_profile_score_and_coverage','{}'))
        rt_gate=bool(profile) and profile['sequence_evalue']<=1e-5 and profile['domain_evalue']<=1e-5 and profile['profile_coverage']>=.6
        result.append(dict(locus_id=lid,dataset='baseline_catalog',family='Type_I',RT_seq=rs[lid],partner_seq=ps[lid],
            pair_complete=r['pair_completeness']=='both_CDS_complete',ART_profile_passes=rt_gate,
            source_DNA_sha256=source_hash(source['path'],source['format'],source['file_sha256']),
            source_path=source['path'],source_format=source['format'],RT_feature_id=cohort[lid]['RT_protein_id'],partner_feature_id=r['partner_id'],
            locus_genbank=str((ROOT/r['result_folder']/'loci'/(lid+'.gb')).relative_to(ROOT)),
            coordinate_map=str((ROOT/r['result_folder']/'loci'/(lid+'.map.json')).relative_to(ROOT)),
            source_genome=r.get('contig_id',lid.split('__')[0]),upstream_nt=int(r['upstream_nt']),
            array_status=r.get('combined_array_assessment','see_baseline'),role=r.get('role','')))
    for folder in [OLD/'deferred_acquisition/processed',OLD/'diverse_partner_search/processed',OUT/'step2/processed',OUT/'step2/parquet/processed']:
        if not (folder/'candidate_RT_proteins.faa').exists():continue
        rs={r.id:str(r.seq) for r in SeqIO.parse(folder/'candidate_RT_proteins.faa','fasta')}
        ps={r.id:str(r.seq) for r in SeqIO.parse(folder/'typeI_partners.faa','fasta')}
        loci={r['locus_id']:r for r in read_tsv(folder/'candidate_loci.tsv')}
        source_meta={v['accession']:v for v in read_tsv(folder/'source_records.tsv')}
        for p in read_tsv(folder/'typeI_pairs.tsv'):
            if p['association']!='primary_nearest_compatible':continue
            r=loci[p['locus_id']]
            result.append(dict(locus_id=p['locus_id'],dataset=str(folder.relative_to(ROOT)),family='Type_I',
                RT_seq=rs[p['RT_key']],partner_seq=ps[p['partner_key']],pair_complete=p['pair_completeness']=='both_CDS_complete',
                ART_profile_passes=bool(json.loads(r['art_profile_score_and_coverage'])),source_genome=p['source_genome'],
                source_DNA_sha256=source_hash(source_meta[p['source_genome']]['path'],source_meta[p['source_genome']]['format'],source_meta[p['source_genome']]['file_sha256']),
                source_path=source_meta[p['source_genome']]['path'],source_format=source_meta[p['source_genome']]['format'],RT_feature_id=r['RT_protein_id'],partner_feature_id=p['partner_id'],
                locus_genbank=str((folder/'loci'/(p['locus_id']+'.gb')).relative_to(ROOT)),coordinate_map=str((folder/'loci'/(p['locus_id']+'.map.json')).relative_to(ROOT)),
                upstream_nt=int(r['upstream_sequence_available']),array_status=r['array_assessment_status'],role=p['role']))
    for name in ['HPGC','JPGD','UHGV_selected']:
        folder=OUT/name
        if not (folder/'RT_proteins.faa').exists():continue
        rs={r.id:str(r.seq) for r in SeqIO.parse(folder/'RT_proteins.faa','fasta')}
        ps={r.id:str(r.seq) for r in SeqIO.parse(folder/'partners.faa','fasta')}
        source_meta={v['accession']:v for v in read_tsv(folder/'acquired_source_records.tsv')}
        for p in read_tsv(folder/'pairs.tsv'):
            if p['association']!='primary_nearest_compatible':continue
            lid=p['locus_id']
            result.append(dict(locus_id=lid,dataset=name,family=p['family'],RT_seq=rs[lid],partner_seq=ps[lid],
                pair_complete=p['pair_complete']=='True',ART_profile_passes=p['ART_profile_passes']=='True',source_genome=p['source_genome'],
                source_DNA_sha256=source_hash(p['source_path'],'genbank',source_meta[p['source_genome']]['file_sha256']),
                source_path=p['source_path'],source_format='genbank',RT_feature_id=p['RT_key'],partner_feature_id=p['partner_key'],locus_genbank=p['locus_genbank'],coordinate_map=p['coordinate_map'],
                upstream_nt=int(p['upstream_nt']),array_status=p['array_status'],role='new_public_catalog',original_source_id=p['original_source_id']))
    # Old type-II additions remain separate from the primary type-I analysis.
    baseline_stage=ROOT/'results/art_expansion/stage1'
    sources={s['accession']:s for s in read_tsv(baseline_stage/'source_records.tsv')}
    for r in read_tsv(baseline_stage/'candidate_loci.tsv'):
        partners=[p for p in json.loads(r['partner_family_and_coordinates']) if p['family']=='Type_II_partner' and p['compatible_adjacent_context']]
        if not partners:continue
        meta=sources[r['source_accession_and_version']]
        rec=SeqIO.read(ROOT/meta['path'].replace('\\','/'),meta['format'])
        feature=lambda pid:next(f for f in rec.features if f.type=='CDS' and f.qualifiers.get('protein_id',[''])[0]==pid)
        rt=feature(r['RT_protein_id']);p=partners[0];pf=feature(p['protein_id'])
        complete=True
        for f in [rt,pf]:
            if re.search(r'[<>]',str(f.location)):complete=False
            try:
                coding=f.extract(rec.seq)[int(f.qualifiers.get('codon_start',['1'])[0])-1:]
                if str(coding.translate(table=int(f.qualifiers.get('transl_table',['11'])[0]),cds=True))!=f.qualifiers['translation'][0]:complete=False
            except Exception:complete=False
        profile=json.loads(r['art_profile_score_and_coverage'])
        result.append(dict(locus_id=r['locus_id'],dataset='baseline_typeII_references',family='Type_II',
            RT_seq=rt.qualifiers['translation'][0],partner_seq=pf.qualifiers['translation'][0],pair_complete=complete,
            ART_profile_passes=profile.get('coverage',0)>=.6,source_genome=rec.id,
            source_DNA_sha256=source_hash(meta['path'],meta['format'],meta['file_sha256']),source_path=meta['path'],source_format=meta['format'],
            RT_feature_id=r['RT_protein_id'],partner_feature_id=p['protein_id'],locus_genbank=r['oriented_locus_path'],coordinate_map=r['coordinate_map_path'],upstream_nt=int(r['upstream_sequence_available']),
            array_status=r['array_assessment_status'],role='baseline_reference'))
    for p in read_tsv(OLD/'panel_screen/new_complete_pair_candidates.tsv'):
        seqs=list(SeqIO.parse(OLD/'panel_screen/new_loci'/(p['locus_id']+'.proteins.faa'),'fasta'))
        meta=next(r for r in read_tsv(OLD/'panel_acquisition.tsv') if r.get('source_accession_version')==p['source_accession_version'])
        result.append(dict(locus_id=p['locus_id'],dataset='prior_typeII_addition',family='Type_II',RT_seq=str(seqs[0].seq),partner_seq=str(seqs[1].seq),
            pair_complete=True,ART_profile_passes=True,source_genome=p['source_accession_version'],source_DNA_sha256=source_hash(meta['path'],meta['format'],meta['file_sha256']),
            source_path=meta['path'],source_format=meta['format'],RT_feature_id=p['RT_protein_id'],partner_feature_id=p['partner_protein_id'],
            locus_genbank=str((OLD/'panel_screen/new_loci'/(p['locus_id']+'.gb')).relative_to(ROOT)),coordinate_map=str((OLD/'panel_screen/new_loci'/(p['locus_id']+'.map.json')).relative_to(ROOT)),
            upstream_nt=int(p['upstream_available_nt']),array_status='not_calibrated',role='new_to_workspace'))
    context={r['run_accession']:r for r in read_tsv(OUT/'step2/parquet/source_run_context.tsv')}
    prior_pairs={hashlib.sha256((r['RT_seq']+'|'+r['partner_seq']).encode()).hexdigest() for r in result if r['dataset'].startswith('baseline') or 'art_sample_expansion_20261104' in r['dataset'] or r['dataset']=='prior_typeII_addition'}
    for r in result:
        matched=re.search(r'(?:SRR|ERR|DRR)\d+',' '.join([r['source_genome'],r['locus_id'],r.get('original_source_id','')]))
        run_id=matched.group() if matched else ''
        info=context.get(run_id,{})
        r.update(source_run=run_id,study_accession=info.get('study_accession','unresolved'),
            source_BioSample=info.get('sample_accession','unresolved'),
            source_group=info.get('sample_accession') or r['source_DNA_sha256'])
        r['protein_pair_present_before_current_run']=hashlib.sha256((r['RT_seq']+'|'+r['partner_seq']).encode()).hexdigest() in prior_pairs
        r['canonical_RT_interval_json']=json.dumps(canonical_RT_interval(r))
    # Exact protein-pair aliases remain in the census but get one analysis vertex.
    # Favor baseline source annotations when catalog reannotation merely changes
    # a partner start on the same DNA and identical RT sequence.
    result.sort(key=lambda r:(not r['dataset'].startswith('baseline'),r['dataset'],r['locus_id']))
    seen={};seen_loci={};seen_regions=defaultdict(list);unique=[];aliases=[]
    for r in result:
        digest=hashlib.sha256((r['RT_seq']+'|'+r['partner_seq']).encode()).hexdigest()
        locus_key=r['source_DNA_sha256']+'|'+hashlib.sha256(r['RT_seq'].encode()).hexdigest()
        if digest in seen:aliases.append(dict(locus_id=r['locus_id'],representative=seen[digest],reason='exact_RT_and_partner_sequence_pair'));continue
        if locus_key in seen_loci:aliases.append(dict(locus_id=r['locus_id'],representative=seen_loci[locus_key],reason='same_source_DNA_and_RT_with_alternative_partner_gene_boundary'));continue
        region=json.loads(r['canonical_RT_interval_json']);equivalent=None
        if region:
            a,b,strand=region
            for other,previous in seen_regions[r['source_DNA_sha256']]:
                x,y,s=previous;overlap=max(0,min(b,y)-max(a,x))
                if strand==s and overlap/min(b-a,y-x)>=.8 and overlap/max(b-a,y-x)>=.8 and ((a-x) if strand==1 else (b-y))%3==0:
                    equivalent=other;break
        if equivalent:
            aliases.append(dict(locus_id=r['locus_id'],representative=equivalent,reason='same_genomic_RT_locus_with_alternative_gene_boundaries'));continue
        seen[digest]=r['locus_id'];r=dict(r,tree_key=f'L{len(unique)+1:05d}',pair_sequence_sha256=digest)
        seen_loci[locus_key]=r['locus_id']
        if region:seen_regions[r['source_DNA_sha256']].append((r['locus_id'],region))
        unique.append(r)
    return result,unique,aliases


def verify_source_maps(rows,folder):
    checks=[]
    for r in rows:
        original=SeqIO.read(ROOT/r['source_path'].replace('\\','/'),r['source_format'])
        locus=SeqIO.read(ROOT/r['locus_genbank'].replace('\\','/'),'genbank')
        mapping=json.loads((ROOT/r['coordinate_map'].replace('\\','/')).read_text())
        pieces=[]
        for s in mapping['segments']:
            text=original.seq[s['source_start0']:s['source_end0']]
            pieces.append(text if s['orientation']==1 else text.reverse_complement())
        assert Seq('').join(pieces)==locus.seq,r['locus_id']
        for idcol,seqcol in [('RT_feature_id','RT_seq'),('partner_feature_id','partner_seq')]:
            f=next(f for f in locus.features if f.type=='CDS' and f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',['']))[0]==r[idcol])
            dna=f.extract(locus.seq)[int(f.qualifiers.get('codon_start',['1'])[0])-1:]
            code=int(f.qualifiers.get('transl_table',['11'])[0])
            if r['pair_complete']:protein=str(dna.translate(table=code,cds=True))
            else:
                protein=str(dna.translate(table=code)).rstrip('*')
                if r[seqcol].startswith('M') and protein[1:]==r[seqcol][1:]:protein='M'+protein[1:]
            assert len(protein)==len(r[seqcol]) and all(a==b or b=='X' for a,b in zip(protein,r[seqcol])),r['locus_id']+' '+idcol
        checks.append(dict(locus_id=r['locus_id'],source_mapping_reconstructed=True,RT_and_partner_translations_replayed=True,pair_complete=r['pair_complete']))
    table(folder/'source_mapping_verification.tsv',checks)
    return dict(pair_source_maps_replayed=len(checks),CDS_translations_replayed=2*len(checks),passed=True)


def tree_analysis(folder,label,records,reference_records=(),reference_labels=None):
    msa=folder/(label+'.aligned.faa');seqpath=folder/(label+'.input.faa')
    records=list(records)+list(reference_records)
    if len(records)<3:
        return dict(status='insufficient_sequences',sequences=len(records)),None
    SeqIO.write(records,seqpath,'fasta')
    run([ROOT/'.venv/bin/mafft','--localpair','--maxiterate','1000','--thread',CPUS,'--inputorder',seqpath],folder/(label+'.mafft.log'),stdout=msa)
    alignment=AlignIO.read(msa,'fasta');columns=[i for i in range(alignment.get_alignment_length()) if sum(r.seq[i]!='-' for r in alignment)>=.5*len(alignment)]
    if len(columns)<50:raise ValueError('Too few alignable columns in '+label)
    trimmed=[SeqRecord(Seq(''.join(str(r.seq[i]) for i in columns)),id=r.id,description='') for r in alignment]
    trim=folder/(label+'.trimmed.faa');SeqIO.write(trimmed,trim,'fasta')
    path=folder/(label+'.FastTree.nwk')
    run([ROOT/'.venv/bin/FastTree','-lg','-gamma','-seed','42',trim],folder/(label+'.FastTree.log'),stdout=path)
    tree=Phylo.read(path,'newick');tips={t.name:t for t in tree.get_terminals()}
    support=[c.confidence for c in tree.get_nonterminals() if c.confidence is not None]
    summary=dict(sequences=len(records),alignment_columns=len(columns),method='MAFFT_L_INS_i_50pct_occupancy_FastTree_LG_Gamma',
        internal_support_kind='FastTree_local_SH_like_support_not_bootstrap',internal_nodes_with_support=len(support),nodes_support_at_least_0_9=sum(v>=.9 for v in support),
        independent_lineages_not_established_by_tree_alone=True)
    return summary,tree


def pair_metrics(a,b,aligner):
    al=aligner.align(a,b)[0]
    first,second=al[0],al[1]
    aligned=sum(x!='-' and y!='-' for x,y in zip(first,second))
    identity=sum(x==y and x!='-' for x,y in zip(first,second))/aligned if aligned else 0
    return identity,min(aligned/len(a),aligned/len(b))


def feasible_blocks(rows,cores,folder):
    aligner=Align.PairwiseAligner(scoring='blastp');aligner.mode='global';aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
    eligible=[];excluded=[]
    for r in rows:
        reasons=[]
        if r['family']!='Type_I':reasons.append('different_partner_family')
        if not r['pair_complete']:reasons.append('incomplete_or_uncertain_CDS')
        if not r['ART_profile_passes']:reasons.append('ART_family_screen_requires_review')
        if r['tree_key'] not in cores:reasons.append('core_boundary_unavailable')
        if reasons:excluded.append(dict(locus_id=r['locus_id'],reasons=';'.join(reasons)))
        else:eligible.append(r)
    edges=defaultdict(set);conflicts=defaultdict(set);comparisons=[]
    for i,a in enumerate(eligible):
        for j in range(i):
            b=eligible[j];ci,cc=pair_metrics(cores[a['tree_key']],cores[b['tree_key']],aligner)
            pi,pc=pair_metrics(a['partner_seq'],b['partner_seq'],aligner)
            ratio=min(len(a['partner_seq']),len(b['partner_seq']))/max(len(a['partner_seq']),len(b['partner_seq']))
            passes=.5<=ci<=.9 and .5<=pi<=.9 and cc>=.8 and pc>=.8 and ratio>=.8 and a['source_group']!=b['source_group']
            if passes:edges[i].add(j);edges[j].add(i)
            if ci>=.95 and pi>=.95 and cc>=.8 and pc>=.8:conflicts[i].add(j);conflicts[j].add(i)
            comparisons.append(dict(locus_a=a['locus_id'],locus_b=b['locus_id'],RT_core_identity=ci,RT_core_mutual_coverage=cc,
                partner_identity=pi,partner_mutual_coverage=pc,partner_shorter_longer_ratio=ratio,passes_sequence_matching=passes))
    triangles=[]
    for i in range(len(eligible)):
        for j in sorted(v for v in edges[i] if v>i):
            for k in sorted(v for v in edges[i]&edges[j] if v>j):triangles.append((i,j,k))
    # Deterministic disjoint packing with a bounded branch-and-bound search.
    # Failure to exhaust the search is reported rather than a maximality claim.
    bynode=defaultdict(list)
    for tri in triangles:
        for v in tri:bynode[v].append(tri)
    best=[];nodes=0;deadline=time.monotonic()+60;timed_out=False
    def search(available,chosen):
        nonlocal best,nodes,timed_out
        nodes+=1
        if time.monotonic()>deadline:timed_out=True;return
        if len(chosen)>len(best):best=list(chosen)
        if len(chosen)+len(available)//3<=len(best):return
        viable={v:[t for t in bynode[v] if all(x in available for x in t)] for v in available}
        options=[v for v in available if viable[v]]
        if not options:return
        v=min(options,key=lambda x:(len(viable[x]),eligible[x]['locus_id']))
        for tri in viable[v]:
            used_sources={eligible[x]['source_group'] for x in tri}
            used_rt={eligible[x]['RT_seq'] for x in tri};used_partner={eligible[x]['partner_seq'] for x in tri}
            near_duplicates=set().union(*(conflicts[x] for x in tri))
            left={x for x in available if x not in tri and x not in near_duplicates and eligible[x]['source_group'] not in used_sources and eligible[x]['RT_seq'] not in used_rt and eligible[x]['partner_seq'] not in used_partner}
            search(left,chosen+[tri])
            if timed_out:return
        search(available-{v},chosen)
    if triangles:search(set(range(len(eligible))),[])
    blocks=[]
    for i,tri in enumerate(best,1):
        blocks.append(dict(block=f'block{i:02d}',loci_json=json.dumps([eligible[v]['locus_id'] for v in tri]),
            status='necessary_sequence_matching_passes_architecture_and_phylogenetic_review_pending',
            source_groups_json=json.dumps([eligible[v]['source_genome'] for v in tri])))
    table(folder/'matching_pair_comparisons.tsv',comparisons);table(folder/'excluded_pairs.tsv',excluded)
    table(folder/'proposed_disjoint_blocks.tsv',blocks)
    table(folder/'all_matching_triples.tsv',[dict(loci_json=json.dumps([eligible[v]['locus_id'] for v in t])) for t in triangles])
    return dict(complete_typeI_pairs_with_core_and_ART_screen=len(eligible),matching_pair_edges=sum(len(v) for v in edges.values())//2,
        matching_triples=len(triangles),disjoint_sequence_matched_blocks_found=len(blocks),packing_search_exhausted=not timed_out,
        packing_search_nodes=nodes,necessary_sequence_gate_at_least_8_blocks=len(blocks)>=8,
        confirmed_eligible_independent_blocks=0,domain_architecture_and_phylogenetic_dependence_review_complete=False,
        confirmatory_structural_design_ready=False,
        interpretation='Sequence-matched candidate blocks are a feasibility upper tier; they do not establish architecture, independence, calibration or power')


def main():
    folder=OUT/'step5';folder.mkdir(exist_ok=True)
    raw,rows,aliases=collect_pairs()
    verification=verify_source_maps(raw,folder)
    write_json(folder/'verification.json',verification)
    table(folder/'exact_pair_aliases.tsv',aliases)
    table(folder/'raw_pair_records.tsv',[{k:v for k,v in r.items() if k not in ['RT_seq','partner_seq']} for r in raw])
    table(folder/'pair_catalog.tsv',[{k:v for k,v in r.items() if k not in ['RT_seq','partner_seq']} for r in rows])
    rts=[SeqRecord(Seq(r['RT_seq']),id=r['tree_key'],description=r['locus_id']) for r in rows]
    SeqIO.write(rts,folder/'RT_proteins.faa','fasta')
    partners=[SeqRecord(Seq(r['partner_seq']),id=r['tree_key'],description=r['locus_id']) for r in rows if r['family']=='Type_I']
    SeqIO.write(partners,folder/'typeI_partners.faa','fasta')
    run([ROOT/'.venv/bin/hmmsearch','--cpu',min(4,CPUS),'--noali','-E','.001','--domE','.001','--domtblout',folder/'RVT_1.domtbl',PROFILES['RVT_1'],folder/'RT_proteins.faa'],folder/'RVT_1.log')
    first={}
    for h in parse_domains(folder/'RVT_1.domtbl'):
        if h['domain_evalue']<=.001 and (h['target'] not in first or h['protein_envelope_start0']<first[h['target']]['protein_envelope_start0']):first[h['target']]=h
    lookup={r['tree_key']:r for r in rows};cores={};bounds=[]
    for key,h in first.items():
        if h['profile_coverage']<.6:continue
        r=lookup[key];a=max(0,h['protein_envelope_start0']-40);b=min(len(r['RT_seq']),h['protein_envelope_end0']+40)
        cores[key]=r['RT_seq'][a:b];bounds.append(dict(tree_key=key,locus_id=r['locus_id'],start0=a,end0=b,definition='first_significant_PF00078_envelope_plus_minus_40aa_operational'))
    table(folder/'core_boundaries.tsv',bounds)
    records=[SeqRecord(Seq(text),id=key,description=lookup[key]['locus_id']) for key,text in cores.items()]
    SeqIO.write(records,folder/'RT_cores.faa','fasta')
    frozen=ROOT/'results/art_family_reconstruction/domain_phylogeny'
    metadata={r['protein_key']:r for r in read_tsv(frozen/'tree_tip_annotations.tsv') if r['role'] in ['protein_only_labelled_reference','published_ART_reference']}
    refs=[r for r in SeqIO.parse(frozen/'domains.faa','fasta') if r.id in metadata]
    rt_summary,rt_tree=tree_analysis(folder,'RT_core',records,refs,metadata)
    partner_summary,partner_tree=tree_analysis(folder,'TypeI_partner',partners)
    placements=[]
    if rt_tree:
        tips={r.name:r for r in rt_tree.get_terminals()}
        for key in cores:
            nearest=sorted((rt_tree.distance(tips[key],tips[ref.id]),ref.id) for ref in refs)[:3]
            labels=[metadata[k]['broad_class'] for d,k in nearest]
            placements.append(dict(tree_key=key,locus_id=lookup[key]['locus_id'],nearest_reference_labels=';'.join(labels),
                nearest_references_json=json.dumps([dict(reference=k,distance=d,label=metadata[k]['broad_class']) for d,k in nearest]),
                placement='ART_reference_affinity_NART_unresolved' if labels.count('ART')>=2 else 'sequence_family_requires_review',
                architecture_used_to_define_placement=False))
    table(folder/'sequence_only_reference_affinity.tsv',placements)
    run([ROOT/'.venv/bin/mmseqs','easy-cluster',folder/'typeI_partners.faa',folder/'partners70',folder/'mmseqs_tmp',
        '--min-seq-id','0.7','-c','0.8','--cov-mode','0','--threads',CPUS,'-k','5','--split-memory-limit','1G'],folder/'mmseqs.log')
    cluster_rows=[line.split('\t') for line in (folder/'partners70_cluster.tsv').read_text().splitlines()]
    groups={a for a,b in cluster_rows}
    cluster_of={member:center for center,member in cluster_rows}
    for r in rows:
        r['partner_cluster70']=cluster_of.get(r['tree_key'],'not_typeI')
        r['sequence_reference_affinity']=next((p['placement'] for p in placements if p['tree_key']==r['tree_key']),'core_unavailable')
        r['new_in_current_run']=r['dataset'] in ['HPGC','JPGD','UHGV_selected'] or 'art_expansion_20261105' in r['dataset']
    exported=lambda seqs:[{k:v for k,v in r.items() if k not in ['RT_seq','partner_seq']} for r in seqs]
    table(folder/'pair_catalog.tsv',exported(rows))
    complete=[r for r in rows if r['family']=='Type_I' and r['pair_complete']]
    strict=[r for r in complete if r['ART_profile_passes'] and r['tree_key'] in cores]
    table(folder/'complete_typeI_pairs.tsv',exported(complete))
    table(folder/'complete_typeI_strict_screen_pairs.tsv',exported(strict))
    representatives=[]
    for cluster in sorted(groups):
        options=[r for r in complete if r['partner_cluster70']==cluster]
        if options:representatives.append(min(options,key=lambda r:(not r['ART_profile_passes'],r['upstream_nt']<6000,-r['upstream_nt'],r['locus_id'])))
    table(folder/'complete_pair_representatives.tsv',exported(representatives))
    for label,chosen in [('complete_typeI',complete),('strict_typeI',strict),('representatives',representatives)]:
        for role,field in [('RT_proteins','RT_seq'),('partners','partner_seq')]:
            SeqIO.write([SeqRecord(Seq(r[field]),id=r['locus_id'],description=r['tree_key']+' provisional_candidate') for r in chosen],folder/(label+'_'+role+'.faa'),'fasta')
    feasibility=feasible_blocks(rows,cores,folder)
    placement_lookup={p['tree_key']:p['placement'] for p in placements}
    sensitivity=[dict(r,ART_profile_passes=r['ART_profile_passes'] or placement_lookup.get(r['tree_key'])=='ART_reference_affinity_NART_unresolved') for r in rows]
    sensdir=folder/'family_review_sensitivity';sensdir.mkdir(exist_ok=True)
    sensitivity_feasibility=feasible_blocks(sensitivity,cores,sensdir)
    summary=dict(raw_pair_records=len(raw),analysis_locus_representatives=len(rows),pair_or_gene_model_aliases=len(aliases),
        typeI_pairs=sum(r['family']=='Type_I' for r in rows),typeII_pairs=sum(r['family']=='Type_II' for r in rows),
        complete_typeI_pairs=len(complete),complete_strict_screen_typeI_pairs=len(strict),complete_partner70_representatives=len(representatives),
        current_run_distinct_typeI_pairs=sum(r['family']=='Type_I' and r['new_in_current_run'] and not r['protein_pair_present_before_current_run'] for r in rows),
        current_run_complete_typeI_pairs=sum(r['new_in_current_run'] and not r['protein_pair_present_before_current_run'] for r in complete),
        typeI_partner70_clusters=len(groups),RT_tree=rt_summary,partner_tree=partner_summary,block_feasibility=feasibility,
        ART_reference_affinity_sensitivity=sensitivity_feasibility,
        sensitivity_limitation='Tree-affinity sensitivity is a review scenario; it does not reclassify low-coverage hits or establish validated ART specificity',
        new_array_significance_tests=False,verified_NART_references=0,verification=verification)
    write_json(folder/'summary.json',summary)
    checkpoint(5,'Separate sequence trees and v2 block-feasibility assessment',summary,
        [folder/'pair_catalog.tsv',folder/'core_boundaries.tsv',folder/'sequence_only_reference_affinity.tsv',folder/'matching_pair_comparisons.tsv',folder/'proposed_disjoint_blocks.tsv',folder/'summary.json'])
    text=['# Five-step ART acquisition run','',
        'All analysis and downloads ran in Slurm CPU batch jobs. Each step has a JSON and Markdown checkpoint in `checkpoints/`.','',
        'The original catalogs and statistical results remain frozen. New rows retain completeness, ART-profile, partner, context, and array evidence as separate fields.','',
        '```json',json.dumps(summary,indent=2),'```','',
        'The trees are descriptive FastTree sequence analyses with local SH-like support, not a completed phylogenetic independence or calibrated specificity study. Proposed blocks pass necessary sequence matching only; architectural review, calibration, simulations and RNA boundaries remain separate gates.','']
    (OUT/'README.md').write_text('\n'.join(text))


if __name__=='__main__':
    import os
    if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use Slurm; login-node analysis is prohibited.')
    main()
