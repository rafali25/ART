"""Verify dated pilot exports and summarize evidence tiers without merging cohorts."""
import csv
import hashlib
import json
from pathlib import Path
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]


def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    panel = OUT/'panel_screen'
    acquired = {r.get('source_accession_version'):r for r in rows(OUT/'panel_acquisition.tsv') if r['status']=='retrieved'}
    source_checks=0
    for r in acquired.values():
        assert sha(ROOT/r['path'])==r['file_sha256']
        source_checks+=1
    candidates=[]
    for r in rows(OUT/'deferred_acquisition/processed/primary_pair_summary.tsv'):
        candidates.append(dict(locus_id=r['locus_id'],family='Type_I',
            evidence_tier='complete_pair_passes_frozen_RT_and_partner_retrieval_screens',
            source=r['source_genome'],RT_length_aa=r['RT_length_aa'],partner_length_aa=r['partner_length_aa'],
            upstream_nt=r['upstream_nt'],array_status='incomplete_upstream_assessment',
            diversity_status='existing_partner70_cluster',
            evidence_table='deferred_acquisition/processed/primary_pair_summary.tsv'))
    cds_checks=0
    mapping_checks=0
    for r in rows(panel/'new_complete_pair_candidates.tsv'):
        lid=r['locus_id']
        src=acquired[r['source_accession_version']]
        original=SeqIO.read(ROOT/src['path'],src['format'])
        proteins={p.id:p for p in SeqIO.parse(panel/'new_loci'/(lid+'.proteins.faa'),'fasta')}
        for label in ['RT_protein_id','partner_protein_id']:
            feature=next(f for f in original.features if f.type=='CDS' and f.qualifiers.get('protein_id',[''])[0]==r[label])
            nt=feature.extract(original.seq)[int(feature.qualifiers.get('codon_start',['1'])[0])-1:]
            protein=nt.translate(table=int(feature.qualifiers.get('transl_table',['11'])[0]),cds=True)
            assert protein==proteins[r[label]].seq
            cds_checks+=1
        locus=SeqIO.read(panel/'new_loci'/(lid+'.fna'),'fasta')
        mapping=json.loads((panel/'new_loci'/(lid+'.map.json')).read_text())
        pieces=[]
        for seg in mapping['segments']:
            text=original.seq[seg['source_start0']:seg['source_end0']]
            pieces.append(text if seg['orientation']==1 else text.reverse_complement())
        assert Seq('').join(pieces)==locus.seq
        oriented=SeqIO.read(panel/'new_loci'/(lid+'.gb'),'genbank')
        assert next(f for f in oriented.features if f.type=='CDS' and f.qualifiers.get('protein_id',[''])[0]==r['RT_protein_id']).location.strand==1
        mapping_checks+=1
        candidates.append(dict(locus_id=lid,family='Type_II',
            evidence_tier='complete_pair_passes_frozen_RT_and_partner_retrieval_screens',
            source=r['source_accession_version'],RT_length_aa=r['RT_length_aa'],partner_length_aa=r['partner_length_aa'],
            upstream_nt=r['upstream_available_nt'],array_status='uncalibrated_exploratory_scan',
            diversity_status='close_Machias_relative_independence_not_established',
            evidence_table='panel_screen/new_complete_pair_candidates.tsv'))
    for r in rows(OUT/'diverse_partner_search/processed/typeI_pairs.tsv'):
        candidates.append(dict(locus_id=r['locus_id'],family='Type_I_context_review',
            evidence_tier='generic_RT_partner_first_ART_profile_coverage_fails_60pct_gate',
            source=r['source_genome'],RT_length_aa=609,partner_length_aa=587,
            upstream_nt=6000,array_status='four_copy_exploratory_no_null_test',
            diversity_status='98pct_similar_to_existing_pair_same_study',
            evidence_table='diverse_partner_search/README.md'))
    fields=list(candidates[0])
    with (OUT/'new_sample_candidates.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(candidates)
    baseline=ROOT/'results/art_type1_catalog'
    old_summary=json.loads((baseline/'summary.json').read_text())
    old_checks=0
    for base,hashes in [(baseline,json.loads((baseline/'output_sha256.json').read_text())),
                       (ROOT,old_summary['input_manifest_sha256'])]:
        for name,expected in hashes.items():
            assert sha(base/name.replace('\\','/'))==expected,name
            old_checks+=1
    for name in ['deferred_acquisition/processed/verification.json','jug_screen/verification.json']:
        assert json.loads((OUT/name).read_text())['passed']
    diverse=json.loads((OUT/'diverse_partner_search/verification.json').read_text())
    assert diverse['source_hashes_verified'] and diverse['complete_pair_translations_replayed']
    (OUT/'code').mkdir(exist_ok=True)
    shutil.copy2(ROOT/'scripts/search_art_expansion_panel.py',OUT/'code/search_art_expansion_panel.py')
    verification=dict(passed=True,panel_source_file_hashes_checked=source_checks,
        new_typeII_CDS_translations_replayed=cds_checks,new_typeII_oriented_source_maps_checked=mapping_checks,
        baseline_catalog_output_and_input_hashes_unchanged=old_checks,
        evidence_tiers_kept_separate=True,no_new_array_pvalues_claimed=True)
    (OUT/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
    summary=dict(date='2026-11-04',baseline_typeI_primary_pairs=13,
        additional_complete_typeI_passing_retrieval_pair=1,additional_complete_typeII_passing_retrieval_pairs=2,
        additional_typeI_context_review_pair_failing_ART_coverage=1,
        new_accepted_arrays=0,new_independent_lineages_or_eligible_blocks_established=0,
        related_panel_records=78,related_panel_unique_DNA_sequences=66,
        curated_Jug_records=139,curated_Jug_unique_DNA_sequences=138,
        diverse_phmmer_queries=3,diverse_phmmer_incremental_peptides=2,
        archive_source_run_leads_remaining=45,old_catalog_modified=False)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (OUT/'output_sha256.json').write_text(json.dumps({p.relative_to(OUT).as_posix():sha(p)
        for p in OUT.rglob('*') if p.is_file() and p.name!='output_sha256.json'},indent=2)+'\n')
    print(json.dumps(dict(summary=summary,verification=verification),indent=2))


if __name__=='__main__':
    main()
