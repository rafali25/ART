"""Audit benchmark coordinates, group exclusions, full-search null replay and p values."""
from collections import defaultdict,Counter
import hashlib
import json

from Bio import SeqIO
from art_array_detector import ROOT, refinement_scan
from art_array_assessment import sensitivity_batch
from extract_art_loci import read_tsv,write_json


def verify():
    out=ROOT/'results/art_validation'
    for path,sha in json.loads((out/'frozen_inputs.json').read_text()).items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha,path
    cfg=json.loads((ROOT/'config/art_array_search.json').read_text())['art_array_search']
    vc=json.loads((ROOT/'config/art_validation.json').read_text())
    rows=read_tsv(out/'sequence_cases.tsv');nulls=read_tsv(out/'null_scores.tsv')
    groups=defaultdict(list)
    for r in nulls:groups[r['case_id'],r['background']].append(r)
    coords=0;replayed=0
    for r in rows:
        rec=SeqIO.read(out/'loci'/(r['case_id']+'.gb'),'genbank')
        rt=next(f for f in rec.features if f.type=='CDS' and f.qualifiers.get('protein_id',[''])[0] in ('synthetic_RT','AAA61471.1','NP_958675.1','ADZ89953.1'))
        start=int(rt.location.start);offset=max(0,start-6000);seq=str(rec.seq[offset:start]).upper()
        assert abs(refinement_scan(seq,cfg,False)['score']-float(r['best_search_score']))<1e-8
        selected=float(r['selected_chain_score']);ps=[]
        for kind,column in (('local_blocks','local_p'),('dinucleotide','dinucleotide_p')):
            draws=groups[r['case_id'],kind]
            assert len(draws)==vc['null_draws_per_background']
            scores=[float(d['score']) for d in draws]
            p=(1+sum(v>=selected for v in scores))/(1+len(scores))
            assert abs(p-float(r[column]))<1e-12
            ps.append(p)
            d=draws[len(draws)//2];seed=int(d['seed'])
            got=sensitivity_batch((seq,kind,cfg,[seed]))[0][1]
            assert abs(got-float(d['score']))<1e-8
            replayed+=1
        assert abs(min(1,max(ps)*len(rows))-float(r['selected_adjusted_p']))<1e-12
        intervals=json.loads(r['detected_intervals_json'])
        assert all(0<=a<b<=start for a,b in intervals)
        assert [c-b for (a,b),(c,d) in zip(intervals,intervals[1:])]==json.loads(r['measured_spacers_json'])
        model=json.loads((out/'details'/(r['case_id']+'.json')).read_text())['model']
        for (a,b),text in zip(intervals,model.get('final_aligned_sequences',[])):
            assert str(rec.seq[a:b]).upper()==text.replace('-','').upper()
            coords+=1
        if r['locus_status']=='supported_ART_locus_candidate':
            assert r['ART_RT_support']==r['partner_support']=='True'
            assert int(r['upstream_shortfall'])==0
            assert r['array_status']=='accepted_array_candidate'
            assert int(r['detected_copies'])>=3
    related={r['protein_id']:r['evaluation_supergroup'] for r in read_tsv(out/'relatedness_groups.tsv')}
    folds=read_tsv(out/'group_holdout_results.tsv')
    for r in folds:
        assert all(related[t]!=r['heldout_group'] for t in r['training_ids'].split(';'))
    rc=read_tsv(out/'reverse_complement_checks.tsv')
    assert all(r[k]=='True' for r in rc for k in ('identical_oriented_DNA','identical_oriented_CDS','equivalent_original_coordinates'))
    result=dict(passed=True,sequence_cases=len(rows),null_scores=len(nulls),independent_null_replays=replayed,
        sequence_coordinate_checks=coords,reverse_complement_loci=len(rc),group_heldout_targets=len(folds),
        behavior_cases_passed=sum(r['pass_expected_behavior']=='True' for r in rows),
        behavior_failures=[r['case_id'] for r in rows if r['pass_expected_behavior']!='True'],
        biological_validation_claimed=False)
    write_json(out/'verification.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':verify()
