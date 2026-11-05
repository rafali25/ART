"""Independent artifact checks for source links, sequence/architecture separation and trees."""
from collections import Counter,defaultdict
import json
from pathlib import Path

from Bio import AlignIO,Phylo,SeqIO
from Bio.Seq import Seq

from art_reconstruction_common import ROOT,sha,table
from art_staged_search import FIELDS
from extract_art_loci import read_tsv,write_json


def verify():
    stage=ROOT/'results/art_expansion/stage1';family=ROOT/'results/art_family_reconstruction'
    sources=read_tsv(stage/'source_records.tsv');records={}
    for r in sources:
        assert sha(ROOT/r['path'])==r['file_sha256'],r['path']
        rec=SeqIO.read(ROOT/r['path'],r['format']);assert rec.id==r['accession'];records[rec.id]=rec
    loci=read_tsv(stage/'candidate_loci.tsv');genes={r['protein_key']:r for r in read_tsv(stage/'protein_to_contig.tsv')}
    architecture={r['protein_key']:r for r in read_tsv(stage/'observed_locus_architecture.tsv')}
    seq={r['protein_key']:r for r in read_tsv(stage/'sequence_family_placement.tsv')}
    placements={r['protein_key']:r for r in read_tsv(family/'domain_phylogeny/sequence_only_placement.tsv')}
    repeats_checked=0;oriented_checked=0
    for r in loci:
        assert set(FIELDS)<=set(r)
        assert r['protein_key'] in genes and r['protein_key'].startswith('g'),'Protein-only reference cannot establish a genomic locus'
        g=genes[r['protein_key']];rec=records[r['contig_id']];f=rec.features[int(g['feature_index'])]
        coords=json.loads(r['rt_coordinates_and_strand'])
        assert (coords['start0'],coords['end0'],coords['strand'])==(int(f.location.start),int(f.location.end),f.location.strand)
        assert r['source_accession_and_version']==rec.id
        oriented=SeqIO.read(ROOT/r['oriented_locus_path'],'genbank');mapping=json.loads((ROOT/r['coordinate_map_path']).read_text())
        pieces=[]
        for s in sorted(mapping['segments'],key=lambda s:s['locus_start0']):
            text=rec.seq[s['source_start0']:s['source_end0']]
            pieces.append(text.reverse_complement() if s['orientation']==-1 else text)
        assert str(Seq('').join(pieces))==str(oriented.seq)
        ort=next(t for t in oriented.features if t.type=='CDS' and t.qualifiers.get('protein_id',[''])[0]==g['protein_id'])
        assert str(ort.extract(oriented.seq))==str(f.extract(rec.seq));oriented_checked+=1
        ntd=json.loads(r['ntd_boundary'])
        if 'end0' in ntd:assert 0<=ntd['end0']<int(g['length_aa'])
        for p in json.loads(r['partner_family_and_coordinates']):
            partner=next(t for t in rec.features if t.type=='CDS' and t.qualifiers.get('protein_id',[''])[0]==p['protein_id'])
            assert int(partner.location.start)==p['start0'] and int(partner.location.end)==p['end0'] and partner.location.strand==p['strand']
        repeats=json.loads(r['repeat_coordinates_and_consensus'])['repeats']
        for rep in repeats:
            text=Seq('').join([rec.seq[s['source_start0']:s['source_end0']].reverse_complement() if s['orientation']==-1 else rec.seq[s['source_start0']:s['source_end0']] for s in rep['source_segments']])
            assert str(text).upper()==rep['sequence']
            assert str(oriented.seq[rep['locus_start0']:rep['locus_end0']]).upper()==rep['sequence'];repeats_checked+=1
        assert [b['locus_start0']-a['locus_end0'] for a,b in zip(repeats,repeats[1:])]==json.loads(r['spacer_lengths'])
        assert seq[r['protein_key']]['placement']==r['phylogenetic_assignment']==placements[r['protein_key']]['assignment']
        assert placements[r['protein_key']]['assignment_uses_architecture']=='False'
        assert 'complete_ART' not in r['final_evidence_category'] and 'supported_ART_locus' not in r['final_evidence_category']
        assert architecture[r['protein_key']]['array_assessment_status']==r['array_assessment_status']
    # Independent union-find replay of DIAMOND's >=140-bit connected component.
    parents={}
    def root(x):
        parents.setdefault(x,x)
        while parents[x]!=x:
            parents[x]=parents[parents[x]];x=parents[x]
        return x
    count=0
    with (family/'replication/all_against_all.tsv').open() as handle:
        for line in handle:
            a,b,bits,*rest=line.split('\t')
            if float(bits)>=140:
                parents[root(b)]=root(a);count+=1
    component=read_tsv(family/'replication/component.tsv');anchor=next(r['protein_key'] for r in component if r['original_id']=='QQM14740.1')
    observed={k for k in parents if root(k)==root(anchor)}
    assert observed=={r['protein_key'] for r in component}
    clusters=defaultdict(list)
    for r in read_tsv(family/'replication/cluster_members.tsv'):clusters[r['mmseqs_cluster']].append(r['member'])
    metadata={r['protein_key']:r for r in read_tsv(family/'protein_metadata.tsv')}
    reps=read_tsv(family/'replication/longest_contig_representatives.tsv')
    for r in reps:
        lengths=[int(metadata[k]['contig_length_nt']) for k in clusters[r['mmseqs_cluster']]]
        assert int(r['contig_length_nt'])==max(lengths)
    fulltree=Phylo.read(family/'replication/full_length.FastTree.nwk','newick')
    assert {t.name for t in fulltree.get_terminals()}=={r['representative'] for r in reps}
    domain=family/'domain_phylogeny';inputs={r.id:r for r in SeqIO.parse(domain/'domain_search_inputs.faa','fasta')}
    fragments={r.id:r for r in SeqIO.parse(domain/'domains.faa','fasta')}
    for r in read_tsv(domain/'domain_coordinates.tsv'):
        assert str(fragments[r['protein_key']].seq)==str(inputs[r['protein_key']].seq[int(r['protein_start0']):int(r['protein_end0'])])
    msa={r.id:str(r.seq) for r in AlignIO.read(domain/'domains.mafft.faa','fasta')}
    columns=json.loads((domain/'alignment_columns.json').read_text())['kept_input_columns0']
    trimmed={r.id:str(r.seq) for r in SeqIO.parse(domain/'domains.occupancy_trimmed.faa','fasta')}
    assert {k:''.join(text[j] for j in columns) for k,text in msa.items()}==trimmed
    tree=Phylo.read(domain/'RT_domain.FastTree.nwk','newick');assert {t.name for t in tree.get_terminals()}==set(trimmed)
    panel=read_tsv(ROOT/'results/art_expansion/stage2/prepared_assembly_sources.tsv')
    for r in panel:
        assert sha(ROOT/r['path'])==r['file_sha256']
        assert sha(ROOT/r['original_annotation_path'])==r['original_annotation_sha256']
        nt=SeqIO.read(ROOT/r['nucleotide_path'],'fasta');annotated=SeqIO.read(ROOT/r['path'],r['format'])
        assert nt.id==annotated.id and str(nt.seq)==str(annotated.seq)
    prep=json.loads((ROOT/'results/art_expansion/stage2/preparation_summary.json').read_text());assert not prep['profile_search_executed']
    assert not json.loads((ROOT/'results/art_expansion/stage3/stage_gate.json').read_text())['allowed']
    result=dict(passed=True,genomic_source_records=len(sources),mapped_candidate_loci=len(loci),oriented_loci_checked=oriented_checked,
        repeat_source_mappings_checked=repeats_checked,DIAMOND_links_replayed=count,sequence_component=len(observed),
        longest_contig_representatives=len(reps),domain_tree_tips=len(trimmed),domain_alignment_columns=len(columns),
        stage2_records_prepared=len(panel),stage2_profile_search_executed=False,stage3_blocked=True,
        sequence_assignment_uses_architecture=False,protein_only_references_exported_as_loci=0)
    write_json(ROOT/'results/art_expansion/verification.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':verify()
