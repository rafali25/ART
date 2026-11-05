"""Separate published-procedure reconstruction and sequence-only RT-domain placement.

The architecture-defined replication clade never supplies the sequence-placement
labels. Protein-only reference tips are not exported as genomic loci.
"""
from __future__ import annotations
import argparse
from collections import defaultdict,Counter
import hashlib
import json
from pathlib import Path
import time

from Bio import AlignIO,Phylo,SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from art_reconstruction_common import ROOT,Runner,sha,table,connected_component,longest_contig_representatives
from extract_art_loci import read_tsv,write_json


def broad_class(label):
    if label.startswith('UG'):return 'UG'
    if label.startswith('Abi'):return 'Abi'
    if label.startswith('G2L'):return 'group_II_like'
    return {'GII':'group_II_intron','Retrons':'retron','DGRs':'DGR','CRISPR':'CRISPR_RT','NotUsed':'unlabelled','UNC':'unclassified'}.get(label,label)


def pool_inputs(stage,out):
    candidates=list(SeqIO.parse(stage/'candidate_RT_proteins.faa','fasta'))
    positions={r['protein_key']:r for r in read_tsv(stage/'sequence_family_placement.tsv')}
    genes={r['protein_key']:r for r in read_tsv(stage/'protein_to_contig.tsv')}
    pool=[];metadata={};excluded=[];outliers=[]
    for r in candidates:
        if not 350<=len(r)<=900:
            excluded.append(dict(protein_key=r.id,reason='outside_published_350_900_replication_range_retained_in_sensitivity_export'));outliers.append(r)
        else:pool.append(r)
        m=positions[r.id]
        metadata[r.id]=dict(protein_key=r.id,original_id=m['RT_protein_id'],broad_class='query',lineage='query',
            role='mapped_candidate',contig_length_nt=int(m['contig_length_nt']),locus_id=m['locus_id'],source_accession=m['source_accession_and_version'],length_aa=len(r))
    for r in SeqIO.parse(ROOT/'.tools/myRT/Models/RTs-collection.faa','fasta'):
        text=str(r.seq).rstrip('*')
        if set(text)-set('ACDEFGHIKLMNPQRSTVWY'):
            excluded.append(dict(protein_key=r.id,reason='noncanonical_reference_protein_not_used'));continue
        key='r'+hashlib.sha256((r.id+'|'+text).encode()).hexdigest()[:16]
        lineage=r.id.rsplit('_',1)[-1]
        pool.append(SeqRecord(Seq(text),id=key,description=r.id))
        metadata[key]=dict(protein_key=key,original_id=r.id,broad_class=broad_class(lineage),lineage=lineage,
            role='protein_only_labelled_reference' if broad_class(lineage) not in ('unlabelled','unclassified') else 'protein_only_unclassified_reference',
            contig_length_nt=0,locus_id='',source_accession='',length_aa=len(text))
    SeqIO.write(pool,out/'pooled_proteins.faa','fasta');SeqIO.write(candidates,out/'all_candidates_including_length_outliers.faa','fasta')
    table(out/'protein_metadata.tsv',list(metadata.values()));table(out/'excluded_or_sensitivity_proteins.tsv',excluded)
    return {r.id:r for r in pool+outliers},metadata,genes


def add_curated_references(pool,metadata,fasta,manifest):
    """Add traceable reference labels, including NART, without locus inference."""
    with Path(fasta).open() as handle:sequences={r.id:r for r in SeqIO.parse(handle,'fasta')}
    added=[]
    for row in read_tsv(manifest):
        if row.get('curation_evidence','').strip()=='' or row.get('source_url','').strip()=='':
            raise ValueError('Curated class references require a source URL and curation evidence')
        r=sequences[row['protein_id']];text=str(r.seq).rstrip('*')
        if set(text)-set('ACDEFGHIKLMNPQRSTVWY'):raise ValueError('Invalid curated reference protein')
        key='v'+hashlib.sha256((r.id+'|'+text).encode()).hexdigest()[:16]
        pool[key]=SeqRecord(Seq(text),id=key,description=r.id)
        metadata[key]=dict(protein_key=key,original_id=r.id,broad_class=row['broad_class'],lineage=row.get('lineage',row['broad_class']),
            role='protein_only_labelled_reference',contig_length_nt=0,locus_id='',source_accession=row.get('source_accession_and_version',''),length_aa=len(text),
            source_url=row['source_url'],curation_evidence=row['curation_evidence'])
        added.append(metadata[key])
    return added


def cluster(run,records,out,prefix,identity,coverage):
    fasta=out/(prefix+'.faa');SeqIO.write(records,fasta,'fasta')
    key=dict(input_sha256=sha(fasta),identity=identity,coverage=coverage,mmseqs_sha256=sha(run.binaries['mmseqs']),kmer_length=5,memory='1G')
    signature=hashlib.sha256(json.dumps(key,sort_keys=True).encode()).hexdigest()[:12]
    stem=out/(prefix+'_'+signature);tmp=out/(prefix+'_'+signature+'_tmp');cache=out/(prefix+'.cache.json')
    result=Path(str(stem)+'_cluster.tsv')
    if cache.exists():
        old=json.loads(cache.read_text())
        if old['key']==key and result.exists() and sha(result)==old['cluster_sha256']:
            return [line.split('\t') for line in result.read_text().splitlines() if line.strip()]
    run.run([run.binaries['mmseqs'],'easy-cluster',fasta,stem,tmp,'--min-seq-id',identity,'-c',coverage,
        '--cov-mode','0','--threads','4','--split-memory-limit','1G','-k','5'],out/(prefix+'.mmseqs.log'),expected=[result])
    write_json(cache,dict(key=key,cluster_sha256=sha(result),cluster_path=str(result)))
    return [line.split('\t') for line in result.read_text().splitlines() if line.strip()]


def replication_workflow(stage,out,pool,metadata,genes,run,members_manifest=None):
    rep=out/'replication';rep.mkdir(exist_ok=True)
    db=rep/'pooled';pairs=rep/'all_against_all.tsv'
    run.run([run.binaries['diamond'],'makedb','--in',out/'pooled_proteins.faa','--db',db],rep/'diamond_makedb.log',expected=[str(db)+'.dmnd'])
    run.run([run.binaries['diamond'],'blastp','--query',out/'pooled_proteins.faa','--db',db,
        '--very-sensitive','--evalue','1e-5','--max-target-seqs','0','--threads','4',
        '--outfmt','6','qseqid','sseqid','bitscore','evalue','qstart','qend','sstart','send','length','pident',
        '--out',pairs],rep/'diamond_blastp.log',expected=[pairs])
    anchor=next((k for k,m in metadata.items() if m['original_id']=='QQM14740.1'),None)
    if anchor is None:raise ValueError('Mapped MarsHill anchor QQM14740.1 is required')
    def links():
        with pairs.open() as handle:
            for line in handle:
                c=line.split('\t')
                if float(c[2])>=140:yield pool[c[0]].id,pool[c[1]].id
    component=connected_component(links(),anchor)
    # The original database's cluster-member lookup is a required external input,
    # not something inferred from our MMseqs 90% clustering below.
    expanded=set(component);extra=[]
    if members_manifest:
        all_records={r.id:r for r in SeqIO.parse(stage/'proteins.faa','fasta')}
        for relation in read_tsv(members_manifest):
            if relation['representative_key'] not in component:continue
            key=relation['member_key']
            if key not in genes or key not in all_records:raise ValueError('Cluster member lacks a retrievable nucleotide mapping: '+key)
            if len(all_records[key])<400:continue
            expanded.add(key);pool[key]=all_records[key]
            if key not in metadata:
                g=genes[key];metadata[key]=dict(protein_key=key,original_id=g['protein_id'],broad_class='query',lineage='query',role='mapped_cluster_member',
                    contig_length_nt=int(g['contig_length_nt']),locus_id=g['accession']+'__'+g['protein_id'],source_accession=g['accession'],length_aa=len(all_records[key]))
            extra.append(dict(representative_key=relation['representative_key'],member_key=key))
    table(rep/'component.tsv',[metadata[k] for k in sorted(component)]);table(rep/'expanded_members.tsv',extra,fields=['representative_key','member_key'])
    SeqIO.write([pool[k] for k in sorted(expanded)],rep/'expanded_component.faa','fasta')
    clusters=cluster(run,[pool[k] for k in sorted(expanded)],rep,'cluster90',.9,.8)
    table(rep/'cluster_members.tsv',[dict(mmseqs_cluster=c,member=m) for c,m in clusters])
    chosen=longest_contig_representatives(clusters,metadata)
    table(rep/'longest_contig_representatives.tsv',[dict(mmseqs_cluster=c,representative=k,contig_length_nt=metadata[k]['contig_length_nt']) for c,k in chosen.items()])
    reps=[pool[k] for k in sorted(chosen.values())];SeqIO.write(reps,rep/'representatives.faa','fasta')
    alignment=rep/'full_length.famsa.faa';treefile=rep/'full_length.FastTree.nwk'
    if len(reps)>=3:
        run.run([run.binaries['famsa'],'-t','4',rep/'representatives.faa',alignment],rep/'famsa.log',expected=[alignment])
        run.run([run.binaries['fasttree'],'-lg','-seed','42',alignment],rep/'FastTree.log',stdout=treefile,expected=[treefile])
    else:write_json(rep/'tree_unavailable.json',dict(reason='fewer_than_three_longest_contig_representatives'))
    arch={r['protein_key']:r for r in read_tsv(stage/'observed_locus_architecture.tsv')}
    anchors=[k for k in chosen.values() if k in arch and arch[k]['architecture_typeI_noncoding500_anchor']=='True']
    architecture_clade=dict(analysis='published_architecture_defined_clade_replication_only',anchors=anchors,
        assignment_not_used_for_sequence_placement=True,rooting_assumption='midpoint rooting of unrooted FastTree, not specified in supplied family methods')
    if len(anchors)>=2 and treefile.exists():
        tree=Phylo.read(treefile,'newick');tree.root_at_midpoint();ancestor=tree.common_ancestor(anchors)
        architecture_clade.update(status='architecture_defined_replication_clade',members=[t.name for t in ancestor.get_terminals()])
        Phylo.write(tree,rep/'full_length.midpoint.nwk','newick')
    else:architecture_clade.update(status='insufficient_independent_typeI_architecture_anchors',members=[])
    write_json(rep/'architecture_defined_clade.json',architecture_clade)
    info=dict(procedure='DIAMOND very-sensitive E<=1e-5; undirected bit-score>=140 MarsHill component; optional database-cluster expansion >=400aa; MMseqs2 90% identity/80% mutual coverage; longest-contig representatives; FAMSA; FastTree LG',
        exact_published_dataset_reproduction=False,component_proteins=len(component),component_reference_proteins=sum(metadata[k]['role'].startswith('protein_only') for k in component),
        expanded_proteins=len(expanded),cluster_representatives=len(reps),architecture_anchors=len(anchors),
        available_baseline_RT_loci=sum(m['role']=='mapped_candidate' and genes[k]['record_role']=='baseline_reference' for k,m in metadata.items() if k in genes),
        published_initial_phage_RT_count=12,published_myRT_reference_count=2019,available_myRT_reference_count=sum(m['role'].startswith('protein_only') for m in metadata.values()),
        original_database_cluster_members_available=bool(members_manifest),original_metagenomic_database_available=False,
        memory_adjustment='MMseqs -k 5 --split-memory-limit 1G on this small Windows panel; the paper does not specify these computational settings',
        input_limitations=['Original 12-seed/2019-reference membership is not supplied','Original 706 added database cluster members and full database representatives are unavailable','Original 823 proteins/230 clusters/95 representatives are not recreated from this bounded panel'])
    write_json(rep/'summary.json',info)
    return info


def reference_panel(stage,out,pool,metadata,run,cfg):
    domain=out/'domain_phylogeny';domain.mkdir(exist_ok=True)
    # All reference families enter selection before any architectural data are read.
    references=[r for k,r in pool.items() if metadata[k]['role']=='protein_only_labelled_reference']
    clusters=cluster(run,references,domain,'reference50',cfg['domain_reference_identity'],cfg['domain_reference_mutual_coverage'])
    table(domain/'reference_cluster_members.tsv',[dict(mmseqs_cluster=c,member=m) for c,m in clusters])
    cluster_of={member:center for center,member in clusters}
    groups=defaultdict(list)
    for r in references:groups[(metadata[r.id]['broad_class'],cluster_of[r.id])].append(r.id)
    representatives=defaultdict(list)
    for (family,clusterid),members in groups.items():
        chosen=sorted(members,key=lambda k:(metadata[k]['lineage'],metadata[k]['original_id']))[0]
        representatives[family].append(chosen)
    selected=[]
    for family,members in sorted(representatives.items()):
        # Deterministic lineage round-robin; no locus features or candidate scores.
        by_lineage=defaultdict(list)
        for k in sorted(members,key=lambda k:metadata[k]['original_id']):by_lineage[metadata[k]['lineage']].append(k)
        ordered=[]
        while any(by_lineage.values()):
            for lineage in sorted(by_lineage):
                if by_lineage[lineage]:ordered.append(by_lineage[lineage].pop(0))
        selected+=ordered[:cfg['domain_reference_limit_per_class']]
    table(domain/'reference_sampling.tsv',[dict(**metadata[k],selection='one_per_50pct_mutual_coverage_cluster_then_lineage_round_robin') for k in selected])
    # Known ART sequence references are duplicated as labelled reference tips;
    # published labels are inherited, not derived from our architecture table.
    art_refs=[]
    for k,m in list(metadata.items()):
        if m['role']=='mapped_candidate' and m['source_accession'] in {r['accession'] for r in read_tsv(ROOT/'results/art_discovery/source_records.tsv')}:
            name='a'+k[1:];art_refs.append(SeqRecord(pool[k].seq,id=name,description='published_ART_reference '+m['original_id']))
            metadata[name]=dict(m,protein_key=name,broad_class='ART',lineage='published_ART',role='published_ART_reference')
    queries=[pool[k] for k,m in metadata.items() if m['role'] in ('mapped_candidate','mapped_cluster_member') and k in pool]
    records=[pool[k] for k in selected]+art_refs+queries
    SeqIO.write(records,domain/'domain_search_inputs.faa','fasta')
    hs,elapsed=run.hmm(run.binaries['RVT_1'],domain/'domain_search_inputs.faa',domain/'RVT_1',z=len(records),report_e=.001)
    table(domain/'RVT_1_hits.tsv',hs)
    first={}
    for h in hs:
        if h['domain_evalue']<=.001 and (h['target'] not in first or h['protein_envelope_start0']<first[h['target']]['protein_envelope_start0']):first[h['target']]=h
    fragments=[];coordinates=[];missing=[]
    for r in records:
        if r.id not in first:missing.append(dict(protein_key=r.id,status='no_significant_RVT_1_domain_not_placed',broad_class=metadata[r.id]['broad_class']));continue
        h=first[r.id];a=max(0,h['protein_envelope_start0']-cfg['domain_flank_aa']);b=min(len(r),h['protein_envelope_end0']+cfg['domain_flank_aa'])
        fragments.append(SeqRecord(r.seq[a:b],id=r.id,description=metadata[r.id]['original_id']))
        coordinates.append(dict(protein_key=r.id,broad_class=metadata[r.id]['broad_class'],role=metadata[r.id]['role'],original_id=metadata[r.id]['original_id'],
            protein_start0=a,protein_end0=b,RVT_1_envelope_start0=h['protein_envelope_start0'],RVT_1_envelope_end0=h['protein_envelope_end0'],
            NTD_start0=0,NTD_end0=h['protein_envelope_start0'],domain_evalue=h['domain_evalue']))
    SeqIO.write(fragments,domain/'domains.faa','fasta');table(domain/'domain_coordinates.tsv',coordinates);table(domain/'unplaced_references.tsv',missing)
    alignment=domain/'domains.mafft.faa'
    run.run(['cmd.exe','/d','/c',ROOT/'.tools/mafft-win/mafft.bat','--localpair','--maxiterate','1000','--thread','4','--inputorder','--out',alignment,domain/'domains.faa'],domain/'MAFFT.log',expected=[alignment])
    msa=AlignIO.read(alignment,'fasta');keep=[j for j in range(msa.get_alignment_length()) if sum(r.seq[j]!='-' for r in msa)>=len(msa)*cfg['domain_alignment_occupancy_trim']]
    trimmed=[SeqRecord(Seq(''.join(str(r.seq[j]) for j in keep)),id=r.id,description='') for r in msa]
    if len(keep)<50:raise ValueError('Too few domain alignment columns to infer a meaningful tree')
    SeqIO.write(trimmed,domain/'domains.occupancy_trimmed.faa','fasta');write_json(domain/'alignment_columns.json',dict(kept_input_columns0=keep,input_columns=msa.get_alignment_length(),occupancy_threshold=cfg['domain_alignment_occupancy_trim']))
    run.run([run.binaries['fasttree'],'-lg','-gamma','-seed','42',domain/'domains.occupancy_trimmed.faa'],domain/'FastTree.log',stdout=domain/'RT_domain.FastTree.nwk',expected=[domain/'RT_domain.FastTree.nwk'])
    tree=Phylo.read(domain/'RT_domain.FastTree.nwk','newick');tips={t.name:t for t in tree.get_terminals()}
    labelled=[k for k in tips if metadata[k]['broad_class'] not in ('query','unlabelled','unclassified')]
    placement=[];nart_count=sum(metadata[k]['broad_class']=='NART' for k in tips)
    for r in queries:
        if r.id not in tips:
            placement.append(dict(protein_key=r.id,assignment='RVT_1_domain_unavailable',nearest_references_json='[]'));continue
        nearest=sorted([(tree.distance(tips[r.id],tips[k]),k) for k in labelled],key=lambda pair:(pair[0],pair[1]))[:5]
        votes=Counter(metadata[k]['broad_class'] for dist,k in nearest[:3]);winner,count=votes.most_common(1)[0]
        assignment='ART_reference_affinity_NART_unresolved' if winner=='ART' and not nart_count else f'{winner}_reference_affinity_provisional'
        placement.append(dict(protein_key=r.id,assignment=assignment,nearest_three_agreement=count/3,
            known_training_or_development_reference=metadata[r.id]['source_accession'] in {r['accession'] for r in read_tsv(ROOT/'results/art_discovery/source_records.tsv')},
            nearest_references_json=json.dumps([dict(protein_key=k,original_id=metadata[k]['original_id'],class_label=metadata[k]['broad_class'],patristic_distance=dist) for dist,k in nearest]),
            NART_comparison_available=bool(nart_count),assignment_uses_architecture=False))
    table(domain/'sequence_only_placement.tsv',placement);table(domain/'tree_tip_annotations.tsv',[metadata[k] for k in tips])
    write_json(domain/'summary.json',dict(method='RVT_1 envelope +/-40aa; MAFFT L-INS-i; 50% occupancy-column trimming; FastTree LG+Gamma approximate phylogeny',
        not_the_published_IQ_TREE_trimAl_domain_analysis=True,columns=len(keep),tips=len(tips),class_counts=dict(Counter(metadata[k]['broad_class'] for k in tips)),
        NART_reference_count=nart_count,scope='Descriptive nearest-labelled-reference affinity, not a calibrated classifier or ART/NART boundary',
        model_limitations=(['No verified NART reference tips'] if not nart_count else [])+['FastTree local support is not IQ-TREE UFBoot/SH-aLRT','Occupancy trim is not trimAl gappyout','The small reference panel does not reproduce the 774-tip study tree']))
    return placement


def merge_sequence_placement(stage,placement):
    lookup={r['protein_key']:r for r in placement}
    seqrows=read_tsv(stage/'sequence_family_placement.tsv')
    for r in seqrows:r['placement']=lookup.get(r['protein_key'],{}).get('assignment','not_in_published_length_reconstruction_range')
    table(stage/'sequence_family_placement.tsv',seqrows)
    rows=read_tsv(stage/'candidate_loci.tsv')
    for r in rows:
        p=lookup.get(r['protein_key'],{})
        r['phylogenetic_assignment']=p.get('assignment','not_placed')
        r['final_evidence_category']='mapped_ART_like_sequence_architecture_provisional' if p.get('assignment','').startswith('ART_reference') else 'mapped_RT_candidate_family_unresolved'
    table(stage/'candidate_loci.tsv',rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',type=Path,default=ROOT/'results/art_expansion/stage1');ap.add_argument('--output',type=Path,default=ROOT/'results/art_family_reconstruction');ap.add_argument('--cluster-members',type=Path);ap.add_argument('--extra-reference-fasta',type=Path);ap.add_argument('--extra-reference-manifest',type=Path);args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);started=time.perf_counter();run=Runner(out)
    cfg=json.loads((args.stage/'run_config.json').read_text());pool,metadata,genes=pool_inputs(args.stage,out)
    if bool(args.extra_reference_fasta)!=bool(args.extra_reference_manifest):raise ValueError('Provide both curated reference FASTA and its evidence manifest')
    if args.extra_reference_fasta:
        added=add_curated_references(pool,metadata,args.extra_reference_fasta,args.extra_reference_manifest)
        table(out/'extra_curated_references.tsv',added)
        table(out/'protein_metadata.tsv',list(metadata.values()))
        SeqIO.write([pool[k] for k,m in metadata.items() if m['role'].startswith('protein_only') or 350<=len(pool[k])<=900],out/'pooled_proteins.faa','fasta')
    if len(pool)>5000:raise ValueError('Bounded reconstruction budget exceeded; measure and configure a larger run separately')
    rep=replication_workflow(args.stage,out,pool,metadata,genes,run,args.cluster_members)
    placement=reference_panel(args.stage,out,pool,metadata,run,cfg);merge_sequence_placement(args.stage,placement)
    write_json(out/'summary.json',dict(replication=rep,domain_query_placements=len(placement),elapsed_seconds=round(time.perf_counter()-started,3),
        sequence_vs_architecture_separated=True,definitive_ART_NART_assignments=0))
    write_json(out/'provenance.json',dict(software=run.metadata,candidate_sequence_sha256=sha(args.stage/'candidate_RT_proteins.faa'),
        myRT_reference_sha256=sha(ROOT/'.tools/myRT/Models/RTs-collection.faa'),config=cfg,
        architecture_used_only_for_separate_replication_clade=True))
    print((out/'summary.json').read_text())


if __name__=='__main__':main()
