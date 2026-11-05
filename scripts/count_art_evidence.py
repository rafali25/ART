"""Read-only evidence census; preserve frozen expansion outputs and decisions."""
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import socket
import statistics

from Bio import SeqIO

from art_expansion_workflow import ROOT, OUT, read_tsv, table, write_json, sha
from art_reconstruction_common import parse_domains


def main():
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run via Slurm; no login-node analysis.')
    target=ROOT/'results/art_evidence_census_20261105'
    target.mkdir(exist_ok=True)
    folder=OUT/'step5'
    rows=read_tsv(folder/'pair_catalog.tsv')
    proteins={r.id:str(r.seq) for r in SeqIO.parse(folder/'RT_proteins.faa','fasta')}
    first={}
    for h in parse_domains(folder/'RVT_1.domtbl'):
        if h['domain_evalue']<=.001 and (h['target'] not in first or h['protein_envelope_start0']<first[h['target']]['protein_envelope_start0']):
            first[h['target']]=h
    baseline={r['locus_id']:r for r in read_tsv(ROOT/'results/art_type1_catalog/typeI_samples.tsv')}
    supported_arrays={r['locus_id'] for r in read_tsv(ROOT/'results/art_type1_catalog/supported_array_candidates.tsv')}
    data=[]
    for row in rows:
        r=dict(row)
        hit=first.get(r['tree_key'])
        domain=bool(hit and hit['profile_coverage']>=.6)
        motif=False;motif_text='';ntd=None
        if domain:
            seq=proteins[r['tree_key']]
            a,b=hit['protein_envelope_start0'],hit['protein_envelope_end0']
            matches=list(re.finditer(r'[YF].DD',seq[a:b]))
            motif=bool(matches)
            motif_text=json.dumps([dict(start0=a+m.start(),sequence=m.group()) for m in matches])
            ntd=a
        locus=SeqIO.read(ROOT/r['locus_genbank'].replace('\\','/'),'genbank')
        def identifier(f):return f.qualifiers.get('protein_id',f.qualifiers.get('locus_tag',['']))[0]
        cds=[f for f in locus.features if f.type=='CDS']
        rt=next(f for f in cds if identifier(f)==r['RT_feature_id'])
        partner=next(f for f in cds if identifier(f)==r['partner_feature_id'])
        assert rt.location.strand==1, r['locus_id']
        a,b=int(rt.location.start),int(rt.location.end)
        others=[f for f in cds if f is not rt]
        overlap=any(int(f.location.start)<a and int(f.location.end)>a for f in others)
        previous=max([int(f.location.end) for f in others if int(f.location.end)<=a],default=0)
        gap=0 if overlap else a-previous
        # This is observed annotation-free DNA, not an estimate from window size.
        observed_500=gap>=500
        pa=int(partner.location.start)
        same_strand=partner.location.strand==1
        intervening=[f for f in cds if f is not rt and f is not partner and
                     int(f.location.start)>=b and int(f.location.start)<pa]
        direct=same_strand and pa>=a and not intervening
        old=baseline.get(r['locus_id'])
        legacy_array=r['locus_id'] in supported_arrays or (r['family']=='Type_II' and r['array_status']=='accepted_array_candidate')
        r.update(ART_profile_screen=r['ART_profile_passes']=='True',
            descriptive_ART_reference_affinity=r['sequence_reference_affinity']=='ART_reference_affinity_NART_unresolved',
            ART_vs_NART_classification_validated=False,
            significant_RVT1_domain=domain,domain_YF_x_DD_pattern=motif,
            RT_domain_and_motif=domain and motif,motif_matches_json=motif_text,
            NTD_before_first_RVT1_aa=ntd,NTD_at_least_100aa_descriptive=ntd is not None and ntd>=100,
            NTD_at_least_150aa_descriptive=ntd is not None and ntd>=150,
            NTD_at_least_180aa_descriptive=ntd is not None and ntd>=180,
            same_strand_partner=same_strand,direct_downstream_partner=direct,
            observed_upstream_annotation_free_nt=gap,
            observed_upstream_noncoding_at_least_500nt=observed_500,
            typeI_direct_partner_and_500nt_anchor=r['family']=='Type_I' and direct and observed_500,
            legacy_supported_array=legacy_array,
            legacy_array_evidence_scope='baseline_typeI_combined_catalog' if r['locus_id'] in supported_arrays else
                'baseline_typeII_reference_assessment' if legacy_array else 'no_supported_legacy_array_recorded',
            array_assessed_in_legacy_typeI_catalog=bool(old),
            array_significance_retested_for_merged_collection=False)
        data.append(r)
    table(target/'sample_evidence.tsv',data)
    flags=['ART_profile_screen','descriptive_ART_reference_affinity','ART_vs_NART_classification_validated',
           'significant_RVT1_domain','domain_YF_x_DD_pattern','RT_domain_and_motif',
           'NTD_at_least_100aa_descriptive','NTD_at_least_150aa_descriptive','NTD_at_least_180aa_descriptive',
           'same_strand_partner','direct_downstream_partner','observed_upstream_noncoding_at_least_500nt',
           'typeI_direct_partner_and_500nt_anchor','legacy_supported_array','array_assessed_in_legacy_typeI_catalog']
    groups={
        'typeI_all':[r for r in data if r['family']=='Type_I'],
        'typeI_complete':[r for r in data if r['family']=='Type_I' and r['pair_complete']=='True'],
        'typeII_all':[r for r in data if r['family']=='Type_II'],
        'typeI_strict_complete':[r for r in data if r['family']=='Type_I' and r['pair_complete']=='True' and r['ART_profile_screen']],
    }
    summaries={}
    for name,chosen in groups.items():
        lengths=[r['NTD_before_first_RVT1_aa'] for r in chosen if r['NTD_before_first_RVT1_aa'] is not None]
        summaries[name]=dict(samples=len(chosen),counts={key:sum(r[key] for r in chosen) for key in flags},
            NTD_length_min=min(lengths) if lengths else None,
            NTD_length_median=statistics.median(lengths) if lengths else None,
            NTD_length_max=max(lengths) if lengths else None,
            array_status_counts=dict(Counter(r['array_status'] for r in chosen)),
            profile_domain_motif_context_overlap=sum(r['ART_profile_screen'] and r['RT_domain_and_motif'] and r['same_strand_partner'] for r in chosen),
            profile_domain_motif_context_legacy_array_overlap=sum(r['ART_profile_screen'] and r['RT_domain_and_motif'] and r['same_strand_partner'] and r['legacy_supported_array'] for r in chosen))
    inputs=[folder/'pair_catalog.tsv',folder/'RT_proteins.faa',folder/'RVT_1.domtbl',ROOT/'results/art_type1_catalog/typeI_samples.tsv',
            ROOT/'results/art_type1_catalog/supported_array_candidates.tsv']
    inputs+=list({ROOT/r['locus_genbank'].replace('\\','/') for r in rows})
    write_json(target/'summary.json',dict(completed_utc=datetime.now(timezone.utc).isoformat(),
        slurm_job_id=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),groups=summaries,
        limitations=[
            'Motif pattern presence inside the RVT_1 envelope is a sequence check, not proof of catalytic activity or correct aligned catalytic position.',
            'NTD thresholds 100/150/180 aa are descriptive; Yoon does not specify a universal long-NTD acceptance cutoff.',
            'No verified NART comparison panel, so descriptive tree affinity is not validated ART-family membership.',
            'Legacy supported arrays retain the old catalog testing-family correction; no merged-cohort significance retest was performed.',
            'The 500nt anchor counts observed annotation-free DNA, not simply the available extraction-window length.',
            'No existing candidate classifications, array calls, or frozen outputs were modified.'],
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},script_sha256=sha(__file__)))
    print(json.dumps(summaries,indent=2),flush=True)


if __name__=='__main__':main()
